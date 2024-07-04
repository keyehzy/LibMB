// Copyright (c) 2024 Matheus Sousa
// SPDX-License-Identifier: BSD-2-Clause

#include "Basis.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <unordered_set>

#include "BosonicBasis.h"
#include "FermionicBasis.h"
#include "GenericBasis.h"
#include "Term.h"

using testing::ElementsAre;
using testing::UnorderedElementsAre;

using enum Operator::Type;
using enum Operator::Statistics;
using enum Operator::Spin;

constexpr int binomial(std::size_t n, std::size_t k) {
  if (k == 0 || k == n) {
    return 1;
  }
  return binomial(n - 1, k - 1) + binomial(n - 1, k);
}

TEST(GenericBasisTest, ConstructorAndAttributes) {
  for (std::size_t orbs = 0; orbs < 6; orbs++) {
    for (std::size_t parts = 0; parts < 6; parts++) {
      GenericBasis basis(orbs, parts);
      EXPECT_EQ(basis.elements().size(), binomial(orbs + parts, parts));
    }
  }
}

TEST(FermionicBasisTest, ConstructorAndAttributes) {
  FermionicBasis basis(2, 2, /*allow_double_occupancy=*/true);

  EXPECT_EQ(basis.particles(), 2);
  EXPECT_EQ(basis.orbitals(), 2);
  EXPECT_EQ(basis.elements().size(), binomial(4, 2));
}

TEST(BasisTest, BasisGeneration) {
  FermionicBasis basis(2, 2, /*allow_double_occupancy=*/true);
  EXPECT_EQ(basis.elements().size(), 6);
  EXPECT_THAT(
      basis.elements_as_vector(), UnorderedElementsAre(
                                      std::vector<Operator>{
                                          Operator::creation<Fermion>(Down, 0),
                                          Operator::creation<Fermion>(Down, 1)},
                                      std::vector<Operator>{
                                          Operator::creation<Fermion>(Up, 1),
                                          Operator::creation<Fermion>(Down, 1)},
                                      std::vector<Operator>{
                                          Operator::creation<Fermion>(Down, 0),
                                          Operator::creation<Fermion>(Up, 1)},
                                      std::vector<Operator>{
                                          Operator::creation<Fermion>(Up, 0),
                                          Operator::creation<Fermion>(Down, 1)},
                                      std::vector<Operator>{
                                          Operator::creation<Fermion>(Up, 0),
                                          Operator::creation<Fermion>(Down, 0)},
                                      std::vector<Operator>{
                                          Operator::creation<Fermion>(Up, 0),
                                          Operator::creation<Fermion>(Up, 1)}));
}

TEST(BasisTest, BosonicBasisGeneration) {
  BosonicBasis basis(2, 2);
  EXPECT_EQ(basis.elements().size(), 3);
  EXPECT_THAT(
      basis.elements_as_vector(), UnorderedElementsAre(
                                      std::vector<Operator>{
                                          Operator::creation<Boson>(Up, 0),
                                          Operator::creation<Boson>(Up, 0)},
                                      std::vector<Operator>{
                                          Operator::creation<Boson>(Up, 0),
                                          Operator::creation<Boson>(Up, 1)},
                                      std::vector<Operator>{
                                          Operator::creation<Boson>(Up, 1),
                                          Operator::creation<Boson>(Up, 1)}));
}

TEST(BasisTest, BasisGenerationDisallowingDoubleOccupation) {
  FermionicBasis basis(2, 1, /*allow_double_occupancy=*/true);
  EXPECT_EQ(basis.elements().size(), 4);
  EXPECT_THAT(
      basis.elements_as_vector(),
      UnorderedElementsAre(
          std::vector<Operator>{Operator::creation<Fermion>(Down, 0)},
          std::vector<Operator>{Operator::creation<Fermion>(Up, 0)},
          std::vector<Operator>{Operator::creation<Fermion>(Down, 1)},
          std::vector<Operator>{Operator::creation<Fermion>(Up, 1)}));
}

TEST(BasisTest, IndexingUnique) {
  FermionicBasis basis(2, 2, /*allow_double_occupancy=*/true);

  std::unordered_set<std::size_t> indices;
  for (const auto& [element, index] : basis.elements()) {
    indices.insert(index);
  }
  EXPECT_EQ(indices.size(), basis.elements().size());

  for (const auto& [element, index] : basis.elements()) {
    EXPECT_GE(index, 0);
    EXPECT_LT(index, basis.elements().size());
  }
}

TEST(BasisTest, IndexingInsideBasis) {
  FermionicBasis basis(2, 2, /*allow_double_occupancy=*/true);

  std::vector<Operator> term = {
      Operator::creation<Fermion>(Up, 0), Operator::creation<Fermion>(Down, 1)};
  EXPECT_TRUE(basis.contains(term));
}

TEST(BasisTest, IndexingOutsideBasis) {
  FermionicBasis basis(2, 2, /*allow_double_occupancy=*/true);

  std::vector<Operator> term = {
      Operator::creation<Fermion>(Up, 0), Operator::creation<Fermion>(Up, 2)};
  EXPECT_FALSE(basis.contains(term));
}

TEST(BasisTest, EqualityOperator) {
  FermionicBasis basis1(2, 2, /*allow_double_occupancy=*/true);
  FermionicBasis basis2(2, 2, /*allow_double_occupancy=*/true);
  FermionicBasis basis3(2, 3, /*allow_double_occupancy=*/true);
  FermionicBasis basis4(3, 2, /*allow_double_occupancy=*/true);

  EXPECT_EQ(basis1, basis2);
  EXPECT_NE(basis1, basis3);
  EXPECT_NE(basis1, basis4);
}

TEST(BasisTest, IndexingEmptyTerm) {
  FermionicBasis basis(2, 2, /*allow_double_occupancy=*/true);

  std::vector<Operator> term = {};
  EXPECT_FALSE(basis.contains(term));
}

TEST(BasisTest, IndexingSingleTermSingleBodyBasis) {
  FermionicBasis basis(1, 1, /*allow_double_occupancy=*/true);

  std::vector<Operator> term = {Operator::creation<Fermion>(Up, 0)};
  EXPECT_TRUE(basis.contains(term));
}

TEST(BasisTest, IndexingSingleTermManyBodyBasis) {
  FermionicBasis basis(2, 2, /*allow_double_occupancy=*/true);

  std::vector<Operator> term = {Operator::creation<Fermion>(Up, 0)};
  EXPECT_FALSE(basis.contains(term));
}

// TEST(DISABLED_BasisTest, SortBasis) {
//   FermionicBasis basis(2, 2, /*allow_double_occupancy=*/true);

//   auto sort_fn = [](const std::vector<Operator>& a,
//                     const std::vector<Operator>& b) {
//     int total_spin_a = 0;
//     for (const auto& op : a) {
//       total_spin_a += static_cast<int>(op.spin());
//     }
//     int total_spin_b = 0;
//     for (const auto& op : b) {
//       total_spin_b += static_cast<int>(op.spin());
//     }
//     return total_spin_a < total_spin_b;
//   };
//   basis.sort(sort_fn);

//   std::vector<Operator> first{
//       Operator::creation<Fermion>(Up, 0), Operator::creation<Fermion>(Up,
//       1)};
//   std::vector<Operator> last{
//       Operator::creation<Fermion>(Down, 0),
//       Operator::creation<Fermion>(Down, 1)};
//   EXPECT_EQ(*basis.elements().begin(), first);
//   EXPECT_EQ(*basis.elements().rbegin(), last);
// }

TEST(PrepareUpAndDownRepresentationTest, EmptyElement) {
  BasisElement element;
  std::vector<int> up(5, 0);
  std::vector<int> down(5, 0);
  prepare_up_and_down_representation(element, up, down);
  EXPECT_EQ(up, std::vector<int>(5, 0));
  EXPECT_EQ(down, std::vector<int>(5, 0));
}

TEST(PrepareUpAndDownRepresentationTest, SingleAnnihilationUp) {
  BasisElement element{Operator::annihilation<Fermion>(Up, 2)};
  std::vector<int> up(5, 0);
  std::vector<int> down(5, 0);
  prepare_up_and_down_representation(element, up, down);
  EXPECT_EQ(up, std::vector<int>({0, 0, -1, 0, 0}));
  EXPECT_EQ(down, std::vector<int>(5, 0));
}

TEST(PrepareUpAndDownRepresentationTest, MultipleOperators) {
  BasisElement element{
      Operator::creation<Fermion>(Down, 1),
      Operator::annihilation<Fermion>(Up, 3),
      Operator::creation<Fermion>(Up, 0)};
  std::vector<int> up(5, 0);
  std::vector<int> down(5, 0);
  prepare_up_and_down_representation(element, up, down);
  EXPECT_EQ(up, std::vector<int>({1, 0, 0, -1, 0}));
  EXPECT_EQ(down, std::vector<int>({0, 1, 0, 0, 0}));
}

TEST(StateStringTest, EmptyElement) {
  BasisElement element;
  FermionicBasis basis(5, 0, /*allow_double_occupancy=*/true);
  std::string expected_state = "|  ,  ,  ,  ,  >";
  std::string actual_state = state_string(element, basis.orbitals());
  EXPECT_EQ(actual_state, expected_state);
}

TEST(StateStringTest, SingleUp) {
  BasisElement element;
  element.push_back(Operator::creation<Fermion>(Up, 1));
  FermionicBasis basis(5, 1, /*allow_double_occupancy=*/true);
  std::string expected_state = "|  ,\u2191 ,  ,  ,  >";
  std::string actual_state = state_string(element, basis.orbitals());
  EXPECT_EQ(actual_state, expected_state);
}

TEST(StateStringTest, MultipleOperators) {
  BasisElement element{
      Operator::creation<Fermion>(Up, 0), Operator::creation<Fermion>(Down, 1),
      Operator::creation<Fermion>(Up, 3)};
  FermionicBasis basis(5, 3, /*allow_double_occupancy=*/true);
  std::string expected_state = "|\u2191 , \u2193,  ,\u2191 ,  >";
  std::string actual_state = state_string(element, basis.orbitals());
  EXPECT_EQ(actual_state, expected_state);
}

TEST(StateStringTest, AllUpAndDown) {
  BasisElement element{
      Operator::creation<Fermion>(Up, 0), Operator::creation<Fermion>(Down, 0),
      Operator::creation<Fermion>(Up, 1), Operator::creation<Fermion>(Down, 1)};
  FermionicBasis basis(5, 4, /*allow_double_occupancy=*/true);
  std::string expected_state = "|\u2191\u2193,\u2191\u2193,  ,  ,  >";
  std::string actual_state = state_string(element, basis.orbitals());
  EXPECT_EQ(actual_state, expected_state);
}
