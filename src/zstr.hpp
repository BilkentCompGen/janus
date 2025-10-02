#pragma once
#include <istream>
#include <ostream>
#include <stdexcept>
#include <vector>
#include <string>
#include <cstring>
#include <zlib.h>

// ---------- input ----------
class gz_streambuf_in : public std::streambuf {
    static constexpr std::size_t kBufSize = 1 << 16; // 64 KiB
    gzFile f_ = nullptr;
    std::vector<char> buf_;
public:
    gz_streambuf_in() : buf_(kBufSize) {}
    ~gz_streambuf_in() override { close(); }

    void open(const char* path) {
        close();
        f_ = gzopen(path, "rb");                   // assume gzip
        if (!f_) throw std::runtime_error("gzopen(rb) failed");
        setg(buf_.data(), buf_.data(), buf_.data());
    }
    void close() {
        if (f_) { gzclose(f_); f_ = nullptr; }
        setg(nullptr, nullptr, nullptr);
    }
    bool is_open() const { return f_ != nullptr; }

protected:
    int_type underflow() override {
        if (!f_) return traits_type::eof();
        int n = gzread(f_, buf_.data(), static_cast<unsigned>(buf_.size()));
        if (n < 0) {
            int err = 0; const char* msg = gzerror(f_, &err);
            throw std::ios_base::failure(msg ? msg : "gzread error");
        }
        if (n == 0) return traits_type::eof();
        setg(buf_.data(), buf_.data(), buf_.data() + n);
        return traits_type::to_int_type(*gptr());
    }
    std::streampos seekoff(std::streamoff, std::ios_base::seekdir, std::ios_base::openmode) override {
        return std::streampos(std::streamoff(-1));
    }
    std::streampos seekpos(std::streampos, std::ios_base::openmode) override {
        return std::streampos(std::streamoff(-1));
    }
};

class gzifstream : public std::istream {
    gz_streambuf_in sbuf_;
public:
    gzifstream() : std::istream(nullptr) { rdbuf(&sbuf_); }
    explicit gzifstream(const char* path) : gzifstream() { open(path); }
    explicit gzifstream(const std::string& path) : gzifstream(path.c_str()) {}

    void open(const char* path) { sbuf_.open(path); clear(); }
    void open(const std::string& path) { open(path.c_str()); }

    bool is_open() const { return sbuf_.is_open(); }
    void close() { sbuf_.close(); clear(); }
};

// ---------- output ----------
class gz_streambuf_out : public std::streambuf {
    static constexpr std::size_t kBufSize = 1 << 16; // 64 KiB
    gzFile f_ = nullptr;
    std::vector<char> buf_;
    int level_ = Z_DEFAULT_COMPRESSION;

    void write_chunk(const char* p, std::size_t n) {
        if (!f_ || n == 0) return;
        int wrote = gzwrite(f_, p, static_cast<unsigned>(n));
        if (wrote <= 0 || static_cast<std::size_t>(wrote) != n) {
            int err = 0; const char* msg = gzerror(f_, &err);
            throw std::ios_base::failure(msg ? msg : "gzwrite error");
        }
    }
public:
    explicit gz_streambuf_out(int level = Z_DEFAULT_COMPRESSION)
        : buf_(kBufSize), level_(level) {
        setp(buf_.data(), buf_.data() + buf_.size());
    }
    ~gz_streambuf_out() override { close(); }

    void open(const char* path) {
        close();
        char mode[8] = "wb";
        if (level_ >= 0 && level_ <= 9) { mode[2] = char('0' + level_); mode[3] = '\0'; }
        f_ = gzopen(path, mode);                  // produce gzip on close
        if (!f_) throw std::runtime_error("gzopen(wb) failed");
        gzbuffer(f_, static_cast<unsigned>(kBufSize)); // zlib internal buffer
        setp(buf_.data(), buf_.data() + buf_.size());
    }

    void set_level(int level) { level_ = level; } // call before open()

    void close() {
        try { sync(); } catch (...) { /* ignore in dtor */ }
        if (f_) { gzclose(f_); f_ = nullptr; }
        setp(nullptr, nullptr);
    }

protected:
    int_type overflow(int_type ch) override {
        if (!f_) return traits_type::eof();
        // flush current buffer
        write_chunk(pbase(), static_cast<std::size_t>(pptr() - pbase()));
        setp(buf_.data(), buf_.data() + buf_.size());
        if (!traits_type::eq_int_type(ch, traits_type::eof())) {
            char c = traits_type::to_char_type(ch);
            write_chunk(&c, 1);
        }
        return traits_type::not_eof(ch);
    }

    std::streamsize xsputn(const char* s, std::streamsize n) override {
        std::streamsize space = epptr() - pptr();
        if (n <= space) {
            std::memcpy(pptr(), s, static_cast<size_t>(n));
            pbump(static_cast<int>(n));
            return n;
        }
        // fill remaining space
        if (space > 0) {
            std::memcpy(pptr(), s, static_cast<size_t>(space));
            pbump(static_cast<int>(space));
        }
        // flush buffer
        overflow(traits_type::eof());
        // write the rest directly
        const char* p = s + space;
        std::streamsize left = n - space;
        write_chunk(p, static_cast<size_t>(left));
        return n;
    }

    int sync() override {
        if (!f_) return 0;
        if (pptr() && pptr() > pbase()) {
            write_chunk(pbase(), static_cast<std::size_t>(pptr() - pbase()));
            setp(buf_.data(), buf_.data() + buf_.size());
        }
        int rc = gzflush(f_, Z_SYNC_FLUSH);
        if (rc != Z_OK && rc != Z_STREAM_END) {
            throw std::ios_base::failure("gzflush error");
        }
        return 0;
    }

    std::streampos seekoff(std::streamoff, std::ios_base::seekdir, std::ios_base::openmode) override {
        return std::streampos(std::streamoff(-1));
    }
    std::streampos seekpos(std::streampos, std::ios_base::openmode) override {
        return std::streampos(std::streamoff(-1));
    }
};

class gzofstream : public std::ostream {
    gz_streambuf_out sbuf_;
public:
    explicit gzofstream(int level = Z_DEFAULT_COMPRESSION)
        : std::ostream(nullptr), sbuf_(level) { rdbuf(&sbuf_); }
    gzofstream(const char* path, int level = Z_DEFAULT_COMPRESSION)
        : gzofstream(level) { open(path); }
    gzofstream(const std::string& path, int level = Z_DEFAULT_COMPRESSION)
        : gzofstream(path.c_str(), level) {}

    void open(const char* path) { sbuf_.open(path); clear(); }
    void open(const std::string& path) { open(path.c_str()); }

    void close() { sbuf_.close(); clear(); }
};
