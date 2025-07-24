/**
 * @file progress_bar.cpp
 * @author Andre Caseiro (andre.v.caseiro@tecnico.ulisboa.pt)
 * @brief Standardized progress bar wrapper implementation
 * @version 0.1
 * @date 2025-07-24
 *
 * @copyright Copyright (c) 2025
 *
 */

#include "progress_bar.h"

#include <indicators/cursor_control.hpp>

ProgressBar::ProgressBar(const std::string& prefix_text, size_t max_progress,
                         size_t update_frequency)
    : _bar{indicators::option::BarWidth{50},
           indicators::option::Start{" ["},
           indicators::option::Fill{"█"},
           indicators::option::Lead{"█"},
           indicators::option::Remainder{"-"},
           indicators::option::End{"]"},
           indicators::option::PrefixText{prefix_text},
           indicators::option::ShowPercentage{true},
           indicators::option::ShowElapsedTime{true},
           indicators::option::ShowRemainingTime{true},
           indicators::option::Stream{std::cerr}},
      _max_progress(max_progress),
      _update_frequency(update_frequency) {
    indicators::show_console_cursor(false);
}

void ProgressBar::update(size_t current) {
    if (current % _update_frequency == 0) {
        force_update(current);
    }
}

void ProgressBar::force_update(size_t current) {
    _bar.set_option(indicators::option::PostfixText{std::to_string(current) + "/" +
                                                    std::to_string(_max_progress)});
    _bar.set_progress(100.0 * current / _max_progress);
}

void ProgressBar::complete() {
    _bar.set_option(indicators::option::PostfixText{std::to_string(_max_progress) + "/" +
                                                    std::to_string(_max_progress)});

    _bar.set_progress(100.0);
    indicators::show_console_cursor(true);
}
