#pragma once

#include <map>
#include <random>  //uniform distribution
#include <set>
#include <string>
#include <vector>

#include "CSV.hpp"
#include "Networks.hpp"
#include "pcg_random.hpp"
#include "stringFormat.hpp"

template <typename T>
struct Node_Epidemic : public Node<T> {
    //* Member variables
    std::string state;
    double transitionRate{0.0};

    //* Generator
    Node_Epidemic() {}
    Node_Epidemic(const T& t_index) : Node<T>(t_index) {}
    Node_Epidemic(const T& t_index, const std::string& t_state) : Node<T>(t_index), state(t_state) {}
};
