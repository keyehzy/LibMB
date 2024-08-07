// Copyright (c) 2024 Matheus Sousa
// SPDX-License-Identifier: BSD-2-Clause

#pragma once

// We encode a single creation/annihilation operator as a byte.
// 0b00000000
//          ^ 0 = creation operator, 1 = annihilation operator
//         ^  0 = boson, 1 = fermion
//        ^   0 = spin up, 1 = spin down
//   ^^^^^      = orbital index (0-31)

#include <bit>
#include <cstdint>
#include <functional>
#include <sstream>
#include <string>
#include <vector>

class Operator {
 public:
  enum class Type { Creation = 0, Annihilation = 1 };
  enum class Statistics { Boson = 0, Fermion = 1 };
  enum class Spin { Up = 0, Down = 1 };

  using UIntType = std::uint8_t;
  static constexpr UIntType Bits = 8 * sizeof(UIntType);
  static constexpr UIntType Operator_Mask = 0x1;                // 0b00000001
  static constexpr UIntType Statistics_Mask = 0x2;              // 0b00000010
  static constexpr UIntType Spin_Mask = 0x4;                    // 0b00000100
  static constexpr UIntType Orbital_Mask = (1 << Bits) - Bits;  // 0b11111000

  constexpr static UIntType max_orbital() {
    return 1 << std::countl_one(Orbital_Mask);
  }

  constexpr Operator(
      Type type, Statistics stats, Spin spin, std::size_t orbital)
      : m_data(static_cast<UIntType>(
            (static_cast<UIntType>(type) << 0) |
            (static_cast<UIntType>(stats) << 1) |
            (static_cast<UIntType>(spin) << 2) |
            (static_cast<UIntType>(orbital) << 3))) {}

  constexpr Operator(const Operator& other) : m_data(other.m_data) {}

  constexpr Operator& operator=(const Operator& other) {
    if (this != &other) {
      m_data = other.m_data;
    }
    return *this;
  }

  constexpr ~Operator() = default;

  constexpr Type type() const {
    return static_cast<Type>(m_data & Operator_Mask);
  }

  constexpr Statistics statistics() const {
    return static_cast<Statistics>((m_data & Statistics_Mask) >> 1);
  }

  constexpr Spin spin() const {
    return static_cast<Spin>((m_data & Spin_Mask) >> 2);
  }

  constexpr std::size_t orbital() const {
    return static_cast<std::size_t>((m_data & Orbital_Mask) >> 3);
  }

  constexpr UIntType identifier() const { return m_data >> 1; }

  constexpr UIntType raw() const { return m_data; }

  constexpr bool operator<(Operator other) const {
    return m_data < other.m_data;
  }

  constexpr bool operator==(Operator other) const {
    return m_data == other.m_data;
  }

  constexpr bool operator!=(Operator other) const {
    return m_data != other.m_data;
  }

  constexpr bool is_boson() const { return statistics() == Statistics::Boson; }

  constexpr bool is_fermion() const {
    return statistics() == Statistics::Fermion;
  }

  friend std::ostream& operator<<(std::ostream& os, Operator op) {
    constexpr auto typeStr = [](Type type) {
      return type == Type::Creation ? "Creation" : "Annihilation";
    };

    constexpr auto spinStr = [](Spin spin) {
      return spin == Spin::Up ? "Up" : "Down";
    };

    os << "Operator { Type: " << typeStr(op.type())
       << ", Spin: " << spinStr(op.spin()) << ", Orbital: " << op.orbital()
       << " }";
    return os;
  }

  constexpr Operator adjoint() const {
    return Operator(
        type() == Type::Creation ? Type::Annihilation : Type::Creation,
        statistics(), spin(), orbital());
  }

  template <Statistics S>
  static constexpr Operator creation(Spin spin, std::size_t orbital) {
    return Operator(Type::Creation, S, spin, orbital);
  }

  template <Statistics S>
  static constexpr Operator annihilation(Spin spin, std::size_t orbital) {
    return Operator(Type::Annihilation, S, spin, orbital);
  }

 private:
  UIntType m_data;
};

template <>
struct std::hash<Operator> {
  size_t operator()(Operator op) const {
    return std::hash<Operator::UIntType>()(op.raw());
  }
};

template <>
struct std::hash<std::vector<Operator>> {
  size_t operator()(const std::vector<Operator>& operators) const {
    size_t hash = 0;
    for (const auto& op : operators) {
      hash_combine(hash, std::hash<Operator>{}(op));
    }
    return hash;
  }

 private:
  template <typename T>
  void hash_combine(size_t& seed, const T& value) const {
    seed ^= std::hash<T>{}(value) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
  }
};
