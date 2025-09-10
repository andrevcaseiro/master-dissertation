#include "time_function.h"

#include <sstream>

std::unique_ptr<TimeFunction> ConstantFunction::clone() const {
    return std::make_unique<ConstantFunction>(*this);
}

std::unique_ptr<TimeFunction> ConstantFunction::operator+(const TimeFunction& other) const {
    return other.addTo(*this);
}

std::unique_ptr<TimeFunction> ConstantFunction::addTo(const ConstantFunction& constant) const {
    return std::make_unique<ConstantFunction>(value + constant.value);
}

std::unique_ptr<TimeFunction> ConstantFunction::addTo(const PulseFunction& pulse) const {
    return std::make_unique<PulseFunction>(pulse.v1 + value, pulse.v2 + value, pulse.td, pulse.tr,
                                           pulse.tf, pulse.pw, pulse.per);
}

std::unique_ptr<TimeFunction> ConstantFunction::addTo(const SumFunction& sum) const {
    auto newSum = std::make_unique<SumFunction>(sum);
    *newSum += value;
    return newSum;
}

TimeFunction& ConstantFunction::operator+=(float v) {
    value += v;
    return *this;
}

TimeFunction& ConstantFunction::operator*=(float v) {
    value *= v;
    return *this;
}

std::string ConstantFunction::to_string() const { return std::to_string(value); }

float PulseFunction::operator()(float t) const {
    /* Initial delay */
    if (t < td) return v1;
    float cycle_time = fmod(t - td, per);

    /* Rising */
    if (cycle_time < tr) return v1 + (v2 - v1) * (cycle_time / tr);

    /* High */
    if (cycle_time < tr + pw) return v2;

    /* Falling */
    if (cycle_time < tr + pw + tf) return v2 - (v2 - v1) * ((cycle_time - tr - pw) / tf);

    /* Low */
    return v1;
}

std::unique_ptr<TimeFunction> PulseFunction::clone() const {
    return std::make_unique<PulseFunction>(*this);
}

std::unique_ptr<TimeFunction> PulseFunction::operator+(const TimeFunction& other) const {
    return other.addTo(*this);
}

std::unique_ptr<TimeFunction> PulseFunction::addTo(const ConstantFunction& constant) const {
    return std::make_unique<PulseFunction>(v1 + constant(0), v2 + constant(0), td, tr, tf, pw, per);
}

std::unique_ptr<TimeFunction> PulseFunction::addTo(const PulseFunction& pulse) const {
    if (td == pulse.td && tr == pulse.tr && tf == pulse.tf && pw == pulse.pw && per == pulse.per) {
        // If the pulse functions have the same time parameters, we can add them directly
        return std::make_unique<PulseFunction>(pulse.v1 + v1, pulse.v2 + v2, pulse.td, pulse.tr,
                                               pulse.tf, pulse.pw, pulse.per);
    } else {
        // If they have different time parameters, create a SumFunction
        auto sumFunc = std::make_unique<SumFunction>();
        sumFunc->addFunction(pulse.clone());
        sumFunc->addFunction(this->clone());
        return sumFunc;
    }
}

std::unique_ptr<TimeFunction> PulseFunction::addTo(const SumFunction& sum) const {
    auto newSum = std::make_unique<SumFunction>(sum);
    newSum->addFunction(this->clone());
    return newSum;
}

TimeFunction& PulseFunction::operator+=(float v) {
    v1 += v;
    v2 += v;
    return *this;
}

TimeFunction& PulseFunction::operator*=(float v) {
    v1 *= v;
    v2 *= v;
    return *this;
}

std::string PulseFunction::to_string() const {
    std::ostringstream oss;
    oss << "pulse(" << v1 << ", " << v2 << ", " << td << ", " << tr << ", " << tf << ", " << pw
        << ", " << per << ")";
    return oss.str();
}

// SumFunction implementation

SumFunction::SumFunction(const SumFunction& other) {
    for (const auto& func : other.functions) {
        functions.push_back(func->clone());
    }
}

SumFunction& SumFunction::operator=(const SumFunction& other) {
    if (this != &other) {
        functions.clear();
        for (const auto& func : other.functions) {
            functions.push_back(func->clone());
        }
    }
    return *this;
}

void SumFunction::addFunction(std::unique_ptr<TimeFunction> func) {
    if (func) {
        functions.push_back(std::move(func));
    }
}

float SumFunction::operator()(float t) const {
    float result = 0.0f;
    for (const auto& func : functions) {
        result += (*func)(t);
    }
    return result;
}

std::unique_ptr<TimeFunction> SumFunction::clone() const {
    return std::make_unique<SumFunction>(*this);
}

std::unique_ptr<TimeFunction> SumFunction::operator+(const TimeFunction& other) const {
    return other.addTo(*this);
}

std::unique_ptr<TimeFunction> SumFunction::addTo(const ConstantFunction& constant) const {
    auto newSum = std::make_unique<SumFunction>(*this);
    *newSum += constant(0);
    return newSum;
}

std::unique_ptr<TimeFunction> SumFunction::addTo(const PulseFunction& pulse) const {
    auto newSum = std::make_unique<SumFunction>(*this);
    newSum->addFunction(pulse.clone());
    return newSum;
}

std::unique_ptr<TimeFunction> SumFunction::addTo(const SumFunction& sum) const {
    auto newSum = std::make_unique<SumFunction>(*this);
    for (const auto& func : sum.functions) {
        newSum->addFunction(func->clone());
    }
    return newSum;
}

TimeFunction& SumFunction::operator+=(float v) {
    if (v != 0.0f) {
        addFunction(std::make_unique<ConstantFunction>(v));
    }
    return *this;
}

TimeFunction& SumFunction::operator*=(float v) {
    for (auto& func : functions) {
        *func *= v;
    }
    return *this;
}

std::string SumFunction::to_string() const {
    std::ostringstream oss;
    oss << "sum(";
    for (size_t i = 0; i < functions.size(); ++i) {
        if (i > 0) oss << " + ";
        oss << functions[i]->to_string();
    }
    oss << ")";
    return oss.str();
}
