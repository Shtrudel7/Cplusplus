#include <iostream>
#include <vector>
#include <algorithm>
#include <limits>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <chrono>
#include <map>
#include <utility>
using namespace std;

struct Item {
    int value;
    int weight;
    Item(int v, int w) : value(v), weight(w) {};
};

bool compareItems(const Item& a, const Item& b) {
    double ratioA = static_cast<double>(a.value) / a.weight;  // Исправлено деление
    double ratioB = static_cast<double>(b.value) / b.weight;
    return ratioA > ratioB;
}


int greedy_KS(int w, vector<Item>& items) {
    sort(items.begin(), items.end(), compareItems);

    int currentWeight = 0;
    int currentValue = 0;
    int maxSingleValue = 0;

    for (const auto& item : items) {
        // Проверяем максимальную ценность одного предмета
        if (item.weight <= w) {
            maxSingleValue = max(maxSingleValue, item.value);
        }

        // Если предмет помещается в рюкзак, добавляем его
        if (currentWeight + item.weight <= w) {
            currentWeight += item.weight;
            currentValue += item.value;
        }
    }
    return max(currentValue, maxSingleValue);
};

int dynamicprog(int capacity, const vector<Item>& items) {
    vector<int> dp(capacity + 1, 0);  // Одномерный массив

    for (const auto& item : items) {
        // Идём справа налево, чтобы не перезаписывать нужные значения
        for (int w = capacity; w >= item.weight; --w) {
            dp[w] = max(dp[w], item.value + dp[w - item.weight]);
        }
    }
    return dp[capacity];
}


double fptas(vector<Item>& items, int w, double epsilon, int n) {
    int F_lb = greedy_KS(w, items);
    double scale = (epsilon * F_lb) / (n * (1 + epsilon));

    vector<Item> scaledItems;
    for (const auto& item : items) {
        int newvalue = static_cast<int>(floor(item.value / scale));
        scaledItems.emplace_back(newvalue, item.weight);
    }

    int approx = dynamicprog(w, scaledItems);
    return floor(approx * scale);
}



int main() {


    int weight, value;
    double epsilon = 0.00001;
    std::vector<Item> items_for200;
    std::ifstream file200("items_100_int2.txt");
    if (!file200.is_open()) {
        std::cerr << "Error opening items_100.txt" << std::endl;
        return 1;
    };

    while (file200 >> value >> weight) {
        items_for200.emplace_back(value, weight);
    };
    file200.close();

    int w200 = 2734;
    int n200 = 100;
  

    std::cout << "=== 100 ITEMS ===" << std::endl;

    // FPTAS for 100 items
    auto start_fptas100 = std::chrono::high_resolution_clock::now();
    double approxvalue = fptas(items_for200, w200, epsilon, n200);
    auto end_fptas100 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_fptas100 = end_fptas100 - start_fptas100;

    std::cout << "FPTAS Algorithm (eps=" << epsilon << "):" << std::endl;
    std::cout << "Approx max value: " << approxvalue << std::endl;
    std::cout << "Time: " << elapsed_fptas100.count() << " seconds" << std::endl;
    std::cout << std::endl;
 
    


    return 0;
}