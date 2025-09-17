// fast_writer_mt.cpp (only the new bits; keep your Segment/BlockInfo/write_node_data)

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <thread>
#include <future>
#include <mutex>
#include <vector>
#include <unistd.h>
#include <cerrno>
#include <system_error>

namespace py = pybind11;

// Robust full-write helper (handles partial pwrite)
static void pwrite_full(int fd, const char* buf, size_t total, off_t pos) {
    size_t done = 0;
    while (done < total) {
        ssize_t n = ::pwrite(fd, buf + done, total - done, pos + done);
        if (n < 0) {
            if (errno == EINTR) continue;
            throw std::system_error(errno, std::generic_category(), "pwrite");
        }
        done += static_cast<size_t>(n);
    }
}

// Same as before but uses pwrite_full
static void write_node_data_full(
    int fd, long long write_pos,
    const std::vector<Segment>& segments,
    int max_read_len, int max_cigar_len)
{
    if (segments.empty()) return;

    size_t record_size = sizeof(int16_t) + max_read_len * 2 + max_cigar_len + sizeof(int16_t) + sizeof(char);
    size_t total_buffer_size = segments.size() * record_size;
    std::vector<char> buffer(total_buffer_size);
    char* cur = buffer.data();

    for (const auto& s : segments) {
        memcpy(cur, &s.offset, sizeof(int16_t)); cur += sizeof(int16_t);

        // seq
        size_t sl = s.seq.size(); if (sl > (size_t)max_read_len) sl = max_read_len;
        memcpy(cur, s.seq.data(), sl);
        memset(cur + sl, 0, max_read_len - sl);
        cur += max_read_len;

        // bq
        size_t bl = s.bq.size(); if (bl > (size_t)max_read_len) bl = max_read_len;
        memcpy(cur, s.bq.data(), bl);
        memset(cur + bl, 0, max_read_len - bl);
        cur += max_read_len;

        // cigar
        size_t cl = s.cigar.size(); if (cl > (size_t)max_cigar_len) cl = max_cigar_len;
        memcpy(cur, s.cigar.data(), cl);
        memset(cur + cl, 0, max_cigar_len - cl);
        cur += max_cigar_len;

        memcpy(cur, &s.rq, sizeof(int16_t)); cur += sizeof(int16_t);
        memcpy(cur, &s.strand, sizeof(char)); cur += sizeof(char);
    }

    pwrite_full(fd, buffer.data(), buffer.size(), write_pos);
}

// Multithreaded flush: parallelize per-node writes
void flush_entire_buffer_mt(
    int fd,
    const std::unordered_map<uint32_t, std::vector<Segment>>& segment_buffer,
    std::unordered_map<uint32_t, BlockInfo>& block_infos,
    const uint32_t block_header_size,
    unsigned max_workers /* e.g., std::thread::hardware_concurrency() */)
{
    if (segment_buffer.empty()) return;
    if (max_workers == 0) max_workers = 1;

    struct Job {
        uint32_t nid;
        const std::vector<Segment>* segs;
        BlockInfo* info;
        long long write_pos;
    };

    std::vector<Job> jobs;
    jobs.reserve(segment_buffer.size());

    // 1) Build jobs with a snapshot of current_pos (serial, no races)
    for (const auto& kv : segment_buffer) {
        uint32_t nid = kv.first;
        auto it = block_infos.find(nid);
        if (it == block_infos.end()) continue;
        BlockInfo* info = &it->second;
        if (kv.second.empty()) continue;

        long long base_offset = static_cast<long long>(info->offset) + block_header_size;
        long long write_pos = base_offset + static_cast<long long>(info->current_pos) * info->record_size;

        jobs.push_back(Job{nid, &kv.second, info, write_pos});
    }

    // 2) Simple thread pool via futures (work stealing could be added if needed)
    std::vector<std::future<void>> futs;
    futs.reserve(jobs.size());
    std::atomic<size_t> idx{0};

    auto worker = [&]() {
        for (;;) {
            size_t i = idx.fetch_add(1);
            if (i >= jobs.size()) break;
            const Job& job = jobs[i];
            write_node_data_full(
                fd, job.write_pos,
                *job.segs,
                job.info->max_read_len,
                job.info->max_cigar_len
            );
        }
    };

    unsigned threads = std::min<unsigned>(max_workers, jobs.size());
    std::vector<std::thread> pool;
    pool.reserve(threads);
    for (unsigned t = 0; t < threads; ++t) pool.emplace_back(worker);
    for (auto& th : pool) th.join();

    // 3) Serially update current_pos after successful writes (avoids races)
    for (const auto& job : jobs) {
        job.info->current_pos += static_cast<uint32_t>(job.segs->size());
    }
}

PYBIND11_MODULE(fast_writer, m) {
    m.doc() = "C++ module for processing and writing segment data buffers.";

    py::class_<Segment>(m, "Segment")
        .def(py::init<int16_t, std::string, std::string, std::string, int16_t, char>());

    py::class_<BlockInfo>(m, "BlockInfo")
        .def(py::init<>())
        .def_readwrite("offset", &BlockInfo::offset)
        .def_readwrite("n_records", &BlockInfo::n_records)
        .def_readwrite("current_pos", &BlockInfo::current_pos)
        .def_readwrite("max_read_len", &BlockInfo::max_read_len)
        .def_readwrite("max_cigar_len", &BlockInfo::max_cigar_len)
        .def_readwrite("record_size", &BlockInfo::record_size)
        .def_readwrite("block_size", &BlockInfo::block_size);

    // Single-thread version (your original)
    m.def("flush_entire_buffer", &flush_entire_buffer,
          py::arg("fd"), py::arg("segment_buffer"),
          py::arg("block_infos"), py::arg("block_header_size"));

    // Multithread version — releases the GIL so threads actually run in parallel
    m.def("flush_entire_buffer_mt",
          &flush_entire_buffer_mt,
          py::arg("fd"), py::arg("segment_buffer"),
          py::arg("block_infos"), py::arg("block_header_size"),
          py::arg("max_workers") = std::thread::hardware_concurrency(),
          py::call_guard<py::gil_scoped_release>());
}
