#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <utility>
#include <set>
#include <tuple>
#include <fstream>
#include <chrono>

// Структура для точки
struct Point {
    double x, y;

    bool operator==(const Point& other) const {
        return ((x == other.x) && (y == other.y));
    };

    bool operator<(const Point& other) const {
        return x < other.x || ((x == other.x) && (y < other.y));
    }

    
};

template <typename T>
T medianOfMedians(const std::vector<T>& nums) {
    // Базовые случаи
    if (nums.empty()) {
        throw std::invalid_argument("Empty");
    }

    switch (nums.size()) {
    case 0:
        throw std::invalid_argument("Array must not be empty");
    case 1:
        return nums[0];
    default:
        if (nums.size() <= 5) {
            // Сортируем массив и возвращаем медиану
            std::vector<T> temp = nums; // Копируем массив
            std::sort(temp.begin(), temp.end());
            return temp[temp.size() / 2];
        }
    }

    // Разбиваем массив на группы по 5 элементов
    std::vector<T> medians;
    for (size_t i = 0; i < nums.size(); i += 5) {
        size_t end = std::min(i + 5, nums.size());
        std::vector<T> chunk(nums.begin() + i, nums.begin() + end);

        // Сортируем группу и берем медиану
        std::sort(chunk.begin(), chunk.end());
        medians.push_back(chunk[chunk.size() / 2]);
    }

    // Рекурсивный вызов для медиан
    return medianOfMedians(medians);
}

// Функция Bridge
std::pair<Point, Point> Bridge(const std::vector<Point>& points, double verticalLine) {
    if (points.size() == 2) {
        Point left = points[0];
        Point right = points[1];
        if (right.x < left.x || (right.x == left.x && right.y < left.y)) {
            std::swap(left, right);
        }
        return {left,right};
    }

    std::vector<Point> candidates;
    std::vector<std::pair<Point, Point>> pairs;
    std::vector<std::pair<std::pair<Point, Point>, double>> slopes; // Храним пары точек и их наклоны

    // Если количество точек нечетное, добавляем последнюю точку в кандидаты
    if (points.size() % 2 != 0) {
        candidates.push_back(points.back());
        //std::cout << "Odd number of points ";
        //points.back().Print();
        
    }

    // Формируем пары
    for (size_t i = 0; i + 1 < points.size(); i += 2) {
        if (points[i].x > points[i + 1].x) {
            pairs.emplace_back(points[i + 1], points[i]); // Меняем местами 
        }
        else {
            pairs.emplace_back(points[i], points[i + 1]);
        }
    }

    // Фильтруем точки с одинаковыми x
    for (auto it = pairs.begin(); it != pairs.end();) {
        if (it->first.x == it->second.x) {
            candidates.push_back(it->first.y > it->second.y ? it->first : it->second);
            it = pairs.erase(it);
        }
        else {
            ++it;
        }
    }

    // Вычисляем углы наклона и добавляем в список с парами
    for (const auto& pair : pairs) {
        double slope = (pair.second.y - pair.first.y) / (pair.second.x - pair.first.x);
        slopes.push_back({ pair, slope });
    }

    // Определяем медианный угол наклона
    std::vector<double> allSlopes;
    for (const auto& slope : slopes) {
        allSlopes.push_back(slope.second);
    };

    double medianSlope = medianOfMedians(allSlopes);

    // Разделяем пары на small, equal, large
    std::vector<std::pair<Point, Point>> small, equal, large;
    for (const auto& slope : slopes) {
        double currentSlope = slope.second;
      

        if (currentSlope < medianSlope) {
            small.push_back(slope.first);
        }
        else if (std::abs(currentSlope - medianSlope) < 1e-9) {
            equal.push_back(slope.first);
        }
        else {
            large.push_back(slope.first);
        }
    }

    // Вычисляем максимальное значение b = y - m*x
    double maxSlope = -1e10;
    Point left, right;
    for (const auto& point : points) {
        double currentSlope = point.y - medianSlope * point.x;
        if (currentSlope > maxSlope) {
            maxSlope = currentSlope;
        }
    }

    // Находим точки left и right
    std::vector<Point> maxSet;
    for (const auto& point : points) {
        if (std::abs(point.y - medianSlope * point.x - maxSlope) < 1e-9) {
            maxSet.push_back(point);
        }
    }
    left = *std::min_element(maxSet.begin(), maxSet.end());
    right = *std::max_element(maxSet.begin(), maxSet.end());

    // Проверяем пересечение с вертикальной линией
    if (left.x <= verticalLine && right.x > verticalLine) {
        return { left, right };
    }

    // Если right.x <= verticalLine
    if (right.x <= verticalLine) {
        for (const auto& pair : large) {
            candidates.push_back(pair.second);
        }
        for (const auto& pair : equal) {
            candidates.push_back(pair.second);
        }
        for (const auto& pair : small) {
            candidates.push_back(pair.first);
            candidates.push_back(pair.second);
        }
    }

    // Если left.x > verticalLine
    if (left.x > verticalLine) {
        for (const auto& pair : small) {
            candidates.push_back(pair.first);
        }
        for (const auto& pair : equal) {
            candidates.push_back(pair.first);
        }
        for (const auto& pair : large) {
            candidates.push_back(pair.first);
            candidates.push_back(pair.second);
        }
    }

    // Рекурсивный вызов для новых кандидатов
    return Bridge(candidates, verticalLine);


};

std::vector<Point> connect(const Point& lower, const Point& upper, const std::vector<Point>& points) {
    if (points.empty()) {
        return {};
    }

    if (lower == upper) {
        return { lower };
    }

    std::vector<double> xCoords;
    for (const auto& point : points) {
        xCoords.push_back(point.x);
    };

    double median = medianOfMedians(xCoords);
    
    auto brd = Bridge(points, median);
    Point left = brd.first;
    Point right = brd.second;
    
    // Разделяем точки на левую и правую части

    std::vector<Point> pointsLeft, pointsRight;

    pointsLeft.push_back(left); 
    pointsRight.push_back(right); 

    for (const auto& point : points) { 
        if (point.x < left.x) { pointsLeft.push_back(point); 
        } 
        else if (point.x > right.x + 1e-9) 
        { pointsRight.push_back(point); } }


    std::vector<Point> output;

    if (left == lower) {
        output.push_back(left);
    }
    else {
        // Рекурсивно обрабатываем левую часть
        auto leftHull = connect(lower, left, pointsLeft);
        output.insert(output.end(), leftHull.begin(), leftHull.end());
    }

    if (right == upper) {
        output.push_back(right);
    }
    else {
        // Рекурсивно обрабатываем правую часть
        auto rightHull = connect(right, upper, pointsRight);
        output.insert(output.end(), rightHull.begin(), rightHull.end());
    };
    // Рекурсивно вычисляем верхнюю оболочку для левой и правой частей
    return output;

};


std::vector<Point> upperHull(const std::vector<Point>& points) {
    if (points.empty()) {
        return {};
    }

    Point minPoint = { std::numeric_limits<double>::max(), std::numeric_limits<double>::min() };
    Point maxPoint = { std::numeric_limits<double>::min(), std::numeric_limits<double>::max() };

    // Находим минимальную и максимальную точки
    for (const auto& point : points) {
        if (point.x < minPoint.x || (point.x == minPoint.x && point.y > minPoint.y)) {
            minPoint = point;
        }
        if (point.x > maxPoint.x || (point.x == maxPoint.x && point.y < maxPoint.y)) {
            maxPoint = point;
        }
    }

    // Если minPoint и maxPoint совпадают, возвращаем одну точку
    if (minPoint == maxPoint) {
        return { minPoint };
    }

    // Формируем временный массив
    std::vector<Point> temporary;
    temporary.push_back(minPoint);
    temporary.push_back(maxPoint);

    for (const auto& point : points) {
        if (point.x > minPoint.x && point.x < maxPoint.x) {
            temporary.push_back(point);
        }
    }

    // Вызываем connect для нахождения верхней оболочки
    return connect(minPoint, maxPoint, temporary);
}

std::vector<Point> lowerHull(std::vector<Point>& points) {
    if (points.empty()) {
        return {};
    }

    // Инвертируем знаки координат "на месте"
    for (auto& point : points) {
        point.x = -point.x;
        point.y = -point.y;
    }

    // Используем ту же логику, что и для верхней оболочки
    auto invertedHull = upperHull(points);

    // Обратно инвертируем координаты "на месте"
    for (auto& point : points) {
        point.x = -point.x;
        point.y = -point.y;
    }

    // Обратно инвертируем координаты в результате upperHull
    for (auto& point : invertedHull) {
        point.x = -point.x;
        point.y = -point.y;
    }

    return invertedHull;
};

std::vector<Point> convexHull(std::vector<Point>& points) {
    if (points.empty()) {
        return {};
    };

    std::vector<Point> upperhull = upperHull(points);
    std::vector<Point> lowerhull = lowerHull(points);
    
    for (const auto& pt : lowerhull) {
        if (std::find(upperhull.begin(), upperhull.end(), pt) == upperhull.end()) {
            upperhull.push_back(pt);
        }
    }

    return upperhull;
}

int main() {

    setlocale(LC_ALL, "RUS");
    std::ifstream inputFile("problem1.txt");
    std::ofstream outputFile("output.txt");

    if (!inputFile.is_open()) {
        std::cerr << "Ошибка: не удалось открыть файл problem.txt" << std::endl;
        return 1;
    }

    if (!outputFile.is_open()) {
        std::cerr << "Ошибка: не удалось открыть файл solution1.txt" << std::endl;
        return 1;
    }

    int n;
    inputFile >> n;

    std::vector<Point> points(n);
    for (int i = 0; i < n; ++i) {
        inputFile >> points[i].x >> points[i].y;
        if (inputFile.fail()) {
            outputFile << "Ошибка ввода данных." << std::endl;
            return 1;
        }
    }
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<Point> result;
    try {
        result = convexHull(points);
    }
    catch (const std::exception& e) {
        outputFile << e.what() << std::endl;
        return 1;
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;


    outputFile << "Количество точек выпуклой оболочки: " << result.size() << std::endl;
    // Вывод результата
    outputFile << "Точки выпуклой оболочки (по часовой стрелке):" << std::endl;
    for (const auto& p : result) {
        outputFile << p.x << " " << p.y << std::endl;
    }
    outputFile << "Время выполнения алгоритма: "
        << elapsed.count() << " секунд." << std::endl;

    inputFile.close();
    outputFile.close();

    return 0;
}


