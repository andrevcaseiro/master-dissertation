/**
 * @file progress_bar.h
 * @author Andre Caseiro (andre.v.caseiro@tecnico.ulisboa.pt)
 * @brief Standardized progress bar wrapper for indicators library
 * @version 0.1
 * @date 2025-07-24
 *
 * @copyright Copyright (c) 2025
 *
 */

#pragma once

#include <indicators/progress_bar.hpp>
#include <string>

/**
 * @brief Standardized progress bar wrapper that provides consistent styling
 * and efficient update handling across the application
 */
class ProgressBar {
   private:
    indicators::ProgressBar _bar;
    size_t _max_progress;
    size_t _update_frequency;

   public:
    /**
     * @brief Construct a new Progress Bar with standardized options
     *
     * @param prefix_text Text to display before the progress bar
     * @param max_progress Maximum value for progress (e.g., number of iterations)
     * @param update_frequency How often to update the display (default: every 100 iterations)
     */
    ProgressBar(const std::string& prefix_text, size_t max_progress, size_t update_frequency = 100);

    /**
     * @brief Update progress with current iteration count
     * Updates display only when needed based on update_frequency
     *
     * @param current Current iteration/progress value
     */
    void update(size_t current);

    /**
     * @brief Force update the progress display regardless of frequency
     *
     * @param current Current iteration/progress value
     */
    void force_update(size_t current);

    /**
     * @brief Mark the progress as completed (100%)
     */
    void complete();

    /**
     * @brief Get the underlying indicators progress bar
     * For advanced usage if needed
     *
     * @return indicators::ProgressBar& Reference to the wrapped progress bar
     */
    indicators::ProgressBar& get_bar() { return _bar; }
};
