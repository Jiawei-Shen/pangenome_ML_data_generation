// fast_writer.cpp
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>          // for BlockInfo in the legacy API
#include <vector>
#include <string>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <unistd.h>                // pwrite (POSIX)
#include <unordered_map>
#include <thread>
#include <mutex>
#include <algorithm>
#include <exception>

namespace py = pybind11;

// ─────────────────────────────────────────────────────────────────────────────
// Data types exposed to Python

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

// Keep the block_infos state in C++ so we don't rebuild it every flush
struct BlockTable {
    std::unordered_map<uint32_t, BlockInfo> m;

    // Build once from Python dict: { nid: {offset, n_records, current_pos, ...} }
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

    // Optional: expose a snapshot of current_pos back to Python if needed later
    py::dict to_py_dict() const {
        py::dict out;
        for (const auto& kv : m) {
            const auto& bi = kv.second;
            py::dict d;
            d["offset"]        = py::int_(bi.offset);
            d["n_records"]     = py::int_(bi.n_records);
            d["current_pos"]   = py::int_(bi.current_pos);
            d["max_read_len"]  = py::int_(bi.max_read_len);
            d["max_cigar_len"] = py::int_(bi.max_cigar_len);
            d["record_size"]   = py::int_(bi.record_size);
            d["block_size"]    = py::int_(bi.block_size);
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

static inline size_t record_size_for(int R, int C) {
    // <h R s R s C s h c>  => 2 + R + R + C + 2 + 1
    return 2u + static_cast<size_t>(R) + static_cast<size_t>(R) + static_cast<size_t>(C) + 2u + 1u;
}

// ─────────────────────────────────────────────────────────────────────────────
// Parallel flush (no per-flush STL conversion; zero-copy Segment access)
//
// Python usage:
//   state = fast_writer.BlockTable(block_infos)  # once
//   fast_writer.flush_entire_buffer_parallel_dict(fd, segment_buffer, state, BLOCK_HDR_SIZE, num_threads)
//
// segment_buffer: dict[int -> list[fast_writer.Segment]]
// state: C++ BlockTable (persistent)

struct Job {
    long long write_pos;
    int R, C;
    std::vector<const Segment*> seg_ptrs;  // raw pointers to C++ Segment (no copies)
    std::vector<py::object> keep_alive;    // holds Python refs so segments don't die
};

void flush_entire_buffer_parallel_dict(
    int fd,
    py::dict segment_buffer,   // Python dict; we iterate it (no big STL build)
    BlockTable& state,
    uint32_t block_header_size,
    int num_threads = (int)std::thread::hardware_concurrency(),
    bool sort_by_offset = true)
{
    // Phase 1 (with GIL): build jobs, compute positions, bump current_pos
    std::vector<Job> jobs;
    jobs.reserve(py::len(segment_buffer));

    for (auto kv : segment_buffer) {
        uint32_t nid = py::cast<uint32_t>(kv.first);
        auto it = state.m.find(nid);
        if (it == state.m.end()) continue;  // unknown nid, skip
        BlockInfo& info = it->second;

        py::list lst = py::cast<py::list>(kv.second);
        const size_t n = py::len(lst);
        if (n == 0) continue;

        // bounds check (optional but recommended)
        if (info.current_pos + n > info.n_records) {
            throw std::runtime_error("Too many records for nid " + std::to_string(nid));
        }

        long long base = static_cast<long long>(info.offset) + block_header_size;
        long long write_pos = base + static_cast<long long>(info.current_pos) * info.record_size;

        Job j;
        j.write_pos = write_pos;
        j.R = static_cast<int>(info.max_read_len);
        j.C = static_cast<int>(info.max_cigar_len);
        j.seg_ptrs.reserve(n);
        j.keep_alive.reserve(n);

        // Collect pointers to underlying C++ Segment without copying
        for (auto obj : lst) {
            py::object o = py::reinterpret_borrow<py::object>(obj);
            // Cast once with GIL held; keep Python ref alive until we finish I/O
            const Segment& s = o.cast<const Segment&>();
            j.keep_alive.emplace_back(std::move(o));
            j.seg_ptrs.push_back(&s);
        }

        // Reserve file space now so next flush computes correct positions
        info.current_pos += static_cast<uint32_t>(n);

        jobs.emplace_back(std::move(j));
    }

    if (jobs.empty()) return;

    if (sort_by_offset) {
        std::sort(jobs.begin(), jobs.end(),
                  [](const Job& a, const Job& b){ return a.write_pos < b.write_pos; });
    }

    if (num_threads <= 0) {
        num_threads = (int)std::thread::hardware_concurrency();
        if (num_threads <= 0) num_threads = 4;
    }
    num_threads = std::min<int>(num_threads, (int)jobs.size());

    std::exception_ptr first_exc = nullptr;
    std::mutex exc_mu;

    // Phase 2 (no GIL): pack & pwrite in parallel
    py::gil_scoped_release nogil;

    auto worker = [&](int tid){
        std::vector<char> scratch; // per-thread buffer
        for (size_t i = tid; i < jobs.size(); i += (size_t)num_threads) {
            const Job& j = jobs[i];

            const size_t rec_sz = record_size_for(j.R, j.C);
            const size_t total  = j.seg_ptrs.size() * rec_sz;
            scratch.resize(total);

            char* p = scratch.data();
            try {
                for (const Segment* s : j.seg_ptrs) {
                    std::memcpy(p, &s->offset, sizeof(int16_t)); p += sizeof(int16_t);
                    cp_field(p, s->seq,   j.R);
                    cp_field(p, s->bq,    j.R);
                    cp_field(p, s->cigar, j.C);
                    std::memcpy(p, &s->rq, sizeof(int16_t)); p += sizeof(int16_t);
                    *p++ = s->strand;
                }

                ssize_t wrote = pwrite(fd, scratch.data(), scratch.size(), j.write_pos);
                if (wrote < 0 || (size_t)wrote != scratch.size()) {
                    throw std::runtime_error("pwrite failed");
                }
            } catch (...) {
                std::lock_guard<std::mutex> lk(exc_mu);
                if (!first_exc) first_exc = std::current_exception();
                // continue to drain; or break;
            }
        }
    };

    std::vector<std::thread> threads;
    threads.reserve(num_threads);
    for (int t = 0; t < num_threads; ++t) threads.emplace_back(worker, t);
    for (auto& th : threads) th.join();

    if (first_exc) std::rethrow_exception(first_exc);
}

// ─────────────────────────────────────────────────────────────────────────────
// Legacy serial API (kept for compatibility).
// NOTE: This version *does* incur pybind11 container conversion on each call.

static void write_node_data(int fd, long long write_pos,
                            const std::vector<Segment>& segments,
                            int R, int C) {
    if (segments.empty()) return;

    const size_t rec_sz = record_size_for(R, C);
    std::vector<char> buffer(segments.size() * rec_sz);
    char* p = buffer.data();

    for (const auto& s : segments) {
        std::memcpy(p, &s.offset, sizeof(int16_t)); p += sizeof(int16_t);
        cp_field(p, s.seq,   R);
        cp_field(p, s.bq,    R);
        cp_field(p, s.cigar, C);
        std::memcpy(p, &s.rq, sizeof(int16_t)); p += sizeof(int16_t);
        *p++ = s.strand;
    }

    ssize_t wrote = pwrite(fd, buffer.data(), buffer.size(), write_pos);
    if (wrote < 0 || (size_t)wrote != buffer.size()) {
        throw std::runtime_error("pwrite failed for a data block.");
    }
}

void flush_entire_buffer( // serial; STL-converting
    int fd,
    const std::unordered_map<uint32_t, std::vector<Segment>>& segment_buffer,
    std::unordered_map<uint32_t, BlockInfo>& block_infos,
    const uint32_t block_header_size)
{
    // (Optional) allow other Python threads; this still pays STL conversion
    py::gil_scoped_release nogil;

    for (const auto& pair : segment_buffer) {
        uint32_t nid = pair.first;
        const std::vector<Segment>& segs = pair.second;
        if (segs.empty()) continue;

        auto it = block_infos.find(nid);
        if (it == block_infos.end()) continue;

        BlockInfo& info = it->second;

        // bounds check
        if (info.current_pos + segs.size() > info.n_records) {
            throw std::runtime_error("too many records for nid " + std::to_string(nid));
        }

        long long base = static_cast<long long>(info.offset) + block_header_size;
        long long write_pos = base + static_cast<long long>(info.current_pos) * info.record_size;

        write_node_data(fd, write_pos, segs, (int)info.max_read_len, (int)info.max_cigar_len);

        info.current_pos += static_cast<uint32_t>(segs.size());
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Module

PYBIND11_MODULE(fast_writer, m) {
    m.doc() = "High-throughput writer for GAM-derived segments with parallel pwrite.";

    py::class_<Segment>(m, "Segment")
        .def(py::init<int16_t, std::string, std::string, std::string, int16_t, char>(),
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

    py::class_(BlockTable)(m, "BlockTable")
        .def(py::init<py::dict>())
        .def("to_py_dict", &BlockTable::to_py_dict);

    // New, fast, parallel, zero-copy path
    m.def("flush_entire_buffer_parallel_dict", &flush_entire_buffer_parallel_dict,
          py::arg("fd"),
          py::arg("segment_buffer"),
          py::arg("state"),
          py::arg("block_header_size"),
          py::arg("num_threads") = (int)std::thread::hardware_concurrency(),
          py::arg("sort_by_offset") = true,
          R"doc(
          Parallel pack & pwrite using Python dict of Segment lists (zero-copy).
          Args:
            fd: open file descriptor for the .dat file
            segment_buffer: dict[int -> list[Segment]]
            state: fast_writer.BlockTable holding block_infos (persistent)
            block_header_size: size in bytes of per-block header (BLOCK_HDR_SIZE)
            num_threads: worker threads (default: hardware_concurrency())
            sort_by_offset: if True, writes in ascending file offset to reduce seeks
          )doc");

    // Legacy serial function (kept for compatibility)
    m.def("flush_entire_buffer", &flush_entire_buffer,
          py::arg("fd"), py::arg("segment_buffer"),
          py::arg("block_infos"), py::arg("block_header_size"),
          R"doc(
          Serial pack & pwrite (STL-converting). Prefer flush_entire_buffer_parallel_dict.
          )doc");
}
