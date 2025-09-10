#pragma once

#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

class ConstantFunction;
class PulseFunction;
class SumFunction;

/**
 * @brief Abstract base class representing a time-dependent function
 *
 */
class TimeFunction {
   public:
    virtual ~TimeFunction() = default;
    virtual float operator()(float t) const = 0;

    // Operator overloads
    virtual std::unique_ptr<TimeFunction> operator+(const TimeFunction& other) const = 0;
    virtual TimeFunction& operator+=(float v) = 0;
    virtual TimeFunction& operator*=(float v) = 0;

    virtual std::unique_ptr<TimeFunction> clone() const = 0;

    virtual std::unique_ptr<TimeFunction> addTo(const ConstantFunction& constant) const = 0;
    virtual std::unique_ptr<TimeFunction> addTo(const PulseFunction& pulse) const = 0;
    virtual std::unique_ptr<TimeFunction> addTo(const SumFunction& sum) const = 0;

    virtual std::string to_string() const = 0;
};

/**
 * @brief Constant function
 *
 */
class ConstantFunction : public TimeFunction {
   private:
    float value;

   public:
    ConstantFunction(float v) : value(v) {}

    float operator()(float /*t*/) const override { return value; }

    std::unique_ptr<TimeFunction> clone() const override;

    std::unique_ptr<TimeFunction> operator+(const TimeFunction& other) const override;
    std::unique_ptr<TimeFunction> addTo(const ConstantFunction& constant) const override;
    std::unique_ptr<TimeFunction> addTo(const PulseFunction& pulse) const override;
    std::unique_ptr<TimeFunction> addTo(const SumFunction& sum) const override;

    TimeFunction& operator+=(float v) override;
    TimeFunction& operator*=(float v) override;

    std::string to_string() const override;
};

/**
 * @brief Pulse function
 *
 */
class PulseFunction : public TimeFunction {
   public:
    float v1, v2, td, tr, tf, pw, per;

    PulseFunction(float v1, float v2, float td, float tr, float tf, float pw, float per)
        : v1(v1), v2(v2), td(td), tr(tr), tf(tf), pw(pw), per(per) {}

    float operator()(float t) const override;

    std::unique_ptr<TimeFunction> clone() const override;

    std::unique_ptr<TimeFunction> operator+(const TimeFunction& other) const override;
    std::unique_ptr<TimeFunction> addTo(const ConstantFunction& constant) const override;
    std::unique_ptr<TimeFunction> addTo(const PulseFunction& pulse) const override;
    std::unique_ptr<TimeFunction> addTo(const SumFunction& sum) const override;

    TimeFunction& operator+=(float v) override;
    TimeFunction& operator*=(float v) override;

    std::string to_string() const override;
};

/**
 * @brief Sum function - represents the sum of multiple time functions
 *
 */
class SumFunction : public TimeFunction {
   private:
    std::vector<std::unique_ptr<TimeFunction>> functions;

   public:
    SumFunction() = default;

    // Constructor that takes a list of functions
    SumFunction(std::vector<std::unique_ptr<TimeFunction>> funcs) : functions(std::move(funcs)) {}

    // Copy constructor
    SumFunction(const SumFunction& other);

    // Assignment operator
    SumFunction& operator=(const SumFunction& other);

    // Add a function to the sum
    void addFunction(std::unique_ptr<TimeFunction> func);

    float operator()(float t) const override;

    std::unique_ptr<TimeFunction> clone() const override;

    std::unique_ptr<TimeFunction> operator+(const TimeFunction& other) const override;
    std::unique_ptr<TimeFunction> addTo(const ConstantFunction& constant) const override;
    std::unique_ptr<TimeFunction> addTo(const PulseFunction& pulse) const override;
    std::unique_ptr<TimeFunction> addTo(const SumFunction& sum) const override;

    TimeFunction& operator+=(float v) override;
    TimeFunction& operator*=(float v) override;

    std::string to_string() const override;
};
