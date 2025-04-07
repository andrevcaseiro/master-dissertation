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
    if (td != pulse.td || tr != pulse.tr || tf != pulse.tf || pw != pulse.pw || per != pulse.per)
        throw std::runtime_error("Can't add pulse functions with different time parameters");

    return std::make_unique<PulseFunction>(pulse.v1 + v1, pulse.v2 + v2, pulse.td, pulse.tr,
                                           pulse.tf, pulse.pw, pulse.per);
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
