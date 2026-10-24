#pragma once
#include <cstdint>
#include <map>
#include <string>
#include <vector>

// The gradient-boosted classifier (bench/train.py, written by bench/export_model.py) and the records of a janus VCF.

struct Node { double value, threshold; int feature, left, right; bool missing_left; };

class Model {
public:
    int features;
    double cut;
    std::map<int, double> cuts;
    explicit Model(const char* path);
    double predict(const double* x) const;
    double threshold(int minimum) const;

private:
    std::vector<double> baseline;
    std::vector<std::vector<std::vector<Node>>> trees;
};

struct Sv {
    std::string line, c1, c2, type;
    int64_t p1, p2;
    uint32_t o;
    std::vector<double> x;
    int own;
    std::vector<int> same, span;
};

class Vcf {
public:
    std::vector<std::string> header;
    std::vector<Sv> r;
    explicit Vcf(const char* path);
    std::vector<std::string> tumors() const;
    int minimum() const;
    void write(const char* path, const std::vector<double>& pr, const std::vector<uint8_t>& du, const std::vector<uint8_t>& pass, const std::vector<std::string>& extra) const;
};

struct Context {
    double tm, nm, c[4];
    Context(const char* hist, int mt);
};

void inputs(const Sv& r, const Context& c, int mt, std::vector<double>& x);
std::vector<uint8_t> duplicates(const std::vector<Sv>& r, const std::vector<double>& pr);
void classify(const Model& m, const std::vector<Sv>& r, const std::vector<size_t>& ix, const char* hist, int mt, std::vector<double>& pr);
