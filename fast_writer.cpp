// fast_writer.cpp
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>       // only for building BlockTable from a Python dict
#include <pybind11/pytypes.h>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>
#include <unordered_map>
#include <stdexcept>
#include <algorithm>
#include <thread>
#include <unistd.h>             // pwrite

namespace py = pybind11;

// ─────────────────────────────────────────────────────────────────────────────
// Data types

struct Segment {
    int16_t offset;
    std::string seq;    // raw bytes
    std::string bq;     // raw bytes
    std::string cigar;  // raw bytes
    int16_t rq;
    char strand;
};

struct BlockInfo {
    uint64_t offset;
    uint32_t n_records;
    uint32_t current_pos;
    uint32_t max_read_len;
    uint32_t max_cigar_len;
    uint32_t record_size;
    uint64_t block_size;
};

// Holds block_infos persistently in C++
struct BlockTable {
    std::unordered_map<uint32_t, BlockInfo> m;

    explicit BlockTable(py::dict py_bi) {
        m.reserve(py::len(py_bi));
        for (auto kv : py_bi) {
            uint32_t nid = py::cast<uint32_t>(kv.first);
            py::dict d   = py::cast<py::dict>(kv.second);

            BlockInfo bi{};
            bi.offset        = py::cast<uint64_t>(d["offset"]);
            bi.n_records     = py::cast<uint32_t>(d["n_records"]);
            bi.current_pos   = py::cast<uint32_t>(d["current_pos"]);
            bi.max_read_len  = py::cast<uint32_t>(d["max_read_len"]);
            bi.max_cigar_len = py::cast<uint32_t>(d["max_cigar_len"]);
            bi.record_size   = py::cast<uint32_t>(d["record_size"]);
            bi.block_size    = py::cast<uint64_t>(d["block_size"]);

            m.emplace(nid, bi);
        }
    }

    py::dict to_py_dict() const {
        py::dict out;
        for (const auto& kv : m) {
            const auto& bi = kv.second;
            py::dict d;
            d["offset"]        = py::int_(static_cast<uint64_t>(bi.offset));
            d["n_records"]     = py::int_(static_cast<uint32_t>(bi.n_records));
            d["current_pos"]   = py::int_(static_cast<uint32_t>(bi.current_pos));
            d["max_read_len"]  = py::int_(static_cast<uint32_t>(bi.max_read_len));
            d["max_cigar_len"] = py::int_(static_cast<uint32_t>(bi.max_cigar_len));
            d["record_size"]   = py::int_(static_cast<uint32_t>(bi.record_size));
            d["block_size"]    = py::int_(static_cast<uint64_t>(bi.block_size));
            out[py::int_(kv.first)] = std::move(d);
        }
        return out;
    }
};

// ─────────────────────────────────────────────────────────────────────────────
// Helpers

static inline void cp_field(char*& p, const std::string& s, int L) {
    size_t n = s.size();
    if (n > static_cast<size_t>(L)) n = static_cast<size_t>(L);
    std::memcpy(p, s.data(), n);
    if (n < static_cast<size_t>(L)) std::memset(p + n, 0, L - n);
    p += L;
}

// <h R s R s C s h c>  => 2 + R + R + C + 2 + 1
static inline size_t record_size_for(int R, int C) {
    return 2u + static_cast<size_t>(R) + static_cast<size_t>(R) + static_cast<size_t>(C) + 2u + 1u;
}

// ─────────────────────────────────────────────────────────────────────────────
// Job struct for parallel pack+write

struct Job {
    uint32_t nid;
    long long write_pos;
    int R, C;
    std::vector<const Segment*> seg_ptrs;  // raw pointers to C++ Segment objects
    std::vector<py::object> keep_alive;    // keep Python wrappers alive
};

// Worker: pack + pwrite (no Python API; GIL is released during whole parallel section)
static void do_job(int fd, const Job& j) {
    const size_t rec_sz = record_size_for(j.R, j.C);
    const size_t total  = j.seg_ptrs.size() * rec_sz;

    std::string buf;
    buf.resize(total);

    char* p = buf.data();
    for (const Segment* s : j.seg_ptrs) {
        std::memcpy(p, &s->offset, sizeof(int16_t)); p += sizeof(int16_t);
        cp_field(p, s->seq,   j.R);
        cp_field(p, s->bq,    j.R);
        cp_field(p, s->cigar, j.C);
        std::memcpy(p, &s->rq, sizeof(int16_t)); p += sizeof(int16_t);
        *p++ = s->strand;
    }

    ssize_t wrote = pwrite(fd, buf.data(), buf.size(), j.write_pos);

    if (wrote < 0 || static_cast<size_t>(wrote) != buf.size()) {
        throw std::runtime_error("pwrite failed for nid " + std::to_string(j.nid));
    }
}

// Main entry: build jobs (with GIL), then release GIL and run parallel pack+writes
void flush_entire_buffer_parallel_dict(
    int fd,
    py::dict segment_buffer,        // dict[int -> list[Segment]]
    BlockTable& state,
    uint32_t block_header_size,
    int num_threads = 4,
    bool sort_by_offset = true
) {
    // 1) Build job list & reserve file space (with GIL)
    std::vector<Job> jobs;
    jobs.reserve(py::len(segment_buffer));

    for (auto kv : segment_buffer) {
        uint32_t nid = py::cast<uint32_t>(kv.first);
        auto it = state.m.find(nid);
        if (it == state.m.end()) continue;  // skip unknown nid

        BlockInfo& info = it->second;
        py::list lst = py::cast<py::list>(kv.second);
        size_t n = py::len(lst);
        if (n == 0) continue;

        if (info.current_pos + n > info.n_records) {
            throw std::runtime_error("Too many records for nid " + std::to_string(nid));
        }

        long long base = static_cast<long long>(info.offset) + block_header_size;
        long long write_pos = base + static_cast<long long>(info.current_pos) * info.record_size;

        Job j;
        j.nid = nid;
        j.write_pos = write_pos;
        j.R = static_cast<int>(info.max_read_len);
        j.C = static_cast<int>(info.max_cigar_len);
        j.seg_ptrs.reserve(n);
        j.keep_alive.reserve(n);

        for (auto obj : lst) {
            py::object o = py::reinterpret_borrow<py::object>(obj);
            const Segment& s = o.cast<const Segment&>();  // C++ ref to bound object
            j.keep_alive.emplace_back(std::move(o));      // keep Python alive until after threads join
            j.seg_ptrs.push_back(&s);
        }

        info.current_pos += static_cast<uint32_t>(n);  // reserve space for next flush
        jobs.emplace_back(std::move(j));
    }

    if (jobs.empty()) return;

    if (sort_by_offset) {
        std::sort(jobs.begin(), jobs.end(),
                  [](const Job& a, const Job& b){ return a.write_pos < b.write_pos; });
    }

    if (num_threads < 1) num_threads = 1;
    if (num_threads > (int)jobs.size()) num_threads = (int)jobs.size();

    // 2) Parallel pack + write with GIL released
    {
        py::gil_scoped_release nogil;

        std::vector<std::thread> workers;
        workers.reserve(num_threads);

        // simple static partition
        auto worker_fn = [&](int tid) {
            size_t start = (jobs.size() * tid) / num_threads;
            size_t end   = (jobs.size() * (tid + 1)) / num_threads;
            for (size_t i = start; i < end; ++i) {
                do_job(fd, jobs[i]);
            }
        };

        for (int t = 0; t < num_threads; ++t) {
            workers.emplace_back(worker_fn, t);
        }
        for (auto& th : workers) th.join();
    }
    // GIL auto-reacquired here; jobs & keep_alive will be destroyed safely.
}

// ─────────────────────────────────────────────────────────────────────────────
// Module

PYBIND11_MODULE(fast_writer, m) {
    m.doc() = "C++ pack+write for GAM-derived segments with parallel pwrite.";

    // Segment: support both str and bytes inputs
    py::class_<Segment>(m, "Segment")
        // str-based ctor
        .def(py::init<int16_t, std::string, std::string, std::string, int16_t, char>(),
             py::arg("offset"), py::arg("seq"), py::arg("bq"),
             py::arg("cigar"), py::arg("rq"), py::arg("strand"))
        // bytes + int strand overload (named args so kwargs work)
        .def(py::init([](int16_t offset,
                         py::bytes seq,
                         py::bytes bq,
                         py::bytes cigar,
                         int16_t rq,
                         int strand) {
                Segment s{};
                s.offset = offset;
                s.seq    = std::string(seq);
                s.bq     = std::string(bq);
                s.cigar  = std::string(cigar);
                s.rq     = rq;
                s.strand = static_cast<char>(strand);
                return s;
            }),
            py::arg("offset"), py::arg("seq"), py::arg("bq"),
            py::arg("cigar"), py::arg("rq"), py::arg("strand"));

    py::class_<BlockInfo>(m, "BlockInfo")
        .def(py::init<>())
        .def_readwrite("offset",        &BlockInfo::offset)
        .def_readwrite("n_records",     &BlockInfo::n_records)
        .def_readwrite("current_pos",   &BlockInfo::current_pos)
        .def_readwrite("max_read_len",  &BlockInfo::max_read_len)
        .def_readwrite("max_cigar_len", &BlockInfo::max_cigar_len)
        .def_readwrite("record_size",   &BlockInfo::record_size)
        .def_readwrite("block_size",    &BlockInfo::block_size);

    py::class_<BlockTable>(m, "BlockTable")
        .def(py::init<py::dict>())
        .def("to_py_dict", &BlockTable::to_py_dict);

    m.def("flush_entire_buffer_parallel_dict", &flush_entire_buffer_parallel_dict,
          py::arg("fd"),
          py::arg("segment_buffer"),
          py::arg("state"),
          py::arg("block_header_size"),
          py::arg("num_threads") = (int)std::thread::hardware_concurrency(),
          py::arg("sort_by_offset") = true,
          R"doc(
              Pack and write all segments in 'segment_buffer' using C++ threads.
              'segment_buffer' is { nid:int -> list[Segment] }.
              'state' is a BlockTable created from your Python block_infos (updates current_pos).
          )doc");
}
