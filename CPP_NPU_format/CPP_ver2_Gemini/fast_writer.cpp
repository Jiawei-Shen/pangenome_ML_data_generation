// file: fast_writer.cpp
#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // For automatic list-to-vector conversion
#include <vector>
#include <string>
#include <cstdint> // For int16_t
#include <cstring> // For memcpy, memset
#include <stdexcept>
#include <unistd.h> // For pwrite (on POSIX systems)

namespace py = pybind11;

// C++ struct to hold segment data from Python
struct Segment {
    int16_t offset;
    std::string seq;
    std::string bq;
    std::string cigar;
    int16_t rq;
    char strand;
};

// The core C++ function to serialize and write data for one node.
// Note: It doesn't need the node ID (nid), but it does need the write position
// and max lengths which are looked up using the nid in the Python script.
void flush_node_buffer(
    int fd,
    long long write_pos,
    const std::vector<Segment>& segments,
    int max_read_len,
    int max_cigar_len)
{
    if (segments.empty()) {
        return;
    }

    // Calculate record and total buffer size
    size_t record_size = sizeof(int16_t) + max_read_len * 2 + max_cigar_len + sizeof(int16_t) + sizeof(char);
    size_t total_buffer_size = segments.size() * record_size;
    std::vector<char> buffer(total_buffer_size);
    char* current_ptr = buffer.data();

    // Serialize all segments into the buffer
    for (const auto& s : segments) {
        memcpy(current_ptr, &s.offset, sizeof(int16_t)); current_ptr += sizeof(int16_t);

        // Copy and pad sequence
        memcpy(current_ptr, s.seq.c_str(), s.seq.length());
        memset(current_ptr + s.seq.length(), 0, max_read_len - s.seq.length());
        current_ptr += max_read_len;

        // Copy and pad base quality
        memcpy(current_ptr, s.bq.c_str(), s.bq.length());
        memset(current_ptr + s.bq.length(), 0, max_read_len - s.bq.length());
        current_ptr += max_read_len;

        // Copy and pad CIGAR
        memcpy(current_ptr, s.cigar.c_str(), s.cigar.length());
        memset(current_ptr + s.cigar.length(), 0, max_cigar_len - s.cigar.length());
        current_ptr += max_cigar_len;

        memcpy(current_ptr, &s.rq, sizeof(int16_t)); current_ptr += sizeof(int16_t);
        memcpy(current_ptr, &s.strand, sizeof(char)); current_ptr += sizeof(char);
    }

    // Write the entire buffer to the file in one operation
    ssize_t bytes_written = pwrite(fd, buffer.data(), buffer.size(), write_pos);
    if (bytes_written == -1) {
        throw std::runtime_error("pwrite failed: " + std::string(strerror(errno)));
    }
    if (static_cast<size_t>(bytes_written) != buffer.size()) {
        throw std::runtime_error("pwrite failed: partial write.");
    }
}

// pybind11 module definition
PYBIND11_MODULE(fast_writer, m) {
    m.doc() = "C++ module for serializing and writing segment data.";

    py::class_<Segment>(m, "Segment")
        // ADD THIS CONSTRUCTOR:
        .def(py::init<int16_t, std::string, std::string, std::string, int16_t, char>(),
            py::arg("offset"), py::arg("seq"), py::arg("bq"),
            py::arg("cigar"), py::arg("rq"), py::arg("strand"))
        .def_readwrite("offset", &Segment::offset)
        .def_readwrite("seq", &Segment::seq)
        .def_readwrite("bq", &Segment::bq)
        .def_readwrite("cigar", &Segment::cigar)
        .def_readwrite("rq", &Segment::rq)
        .def_readwrite("strand", &Segment::strand);

    m.def("flush_node_buffer", &flush_node_buffer,
          "Serializes and writes a list of Segments to a file.",
          py::arg("fd"), py::arg("write_pos"), py::arg("segments"),
          py::arg("max_read_len"), py::arg("max_cigar_len"));
}