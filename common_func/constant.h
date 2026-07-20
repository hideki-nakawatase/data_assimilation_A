#ifndef CONSTANT_H
#define CONSTANT_H
#pragma once

inline constexpr int days = 730;
inline constexpr double dt = 0.005;
inline constexpr int time_steps = 0.2 / dt * days;
inline constexpr int N = 40;
inline constexpr int M = 40;
inline constexpr double F = 8;
inline constexpr double delta = 1.05;
inline constexpr double ϵ = 1e-5;
inline constexpr double alpha = 0.27;

#endif