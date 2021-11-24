#include <iostream>
#include <cmath>
#include <algorithm>
#include <fstream>

constexpr double ERR = 1e-5;

struct Point
{
    double x_;
    double y_;

    Point(double x, double y) : x_{x}, y_{y} {}
};

double distance(const Point& p1, const Point& p2)
{
    return sqrt((p1.x_ - p2.x_) * (p1.x_ - p2.x_) + (p1.y_ - p2.y_) * (p1.y_ - p2.y_));
}

Point next_point_MPI(Point p)
{
    return Point{cos(p.y_) + 0.85, sin(p.x_) - 1.32};
}

Point f_Newton(const Point& p)
{
    return Point{cos(p.y_) - p.x_ + 0.85, sin(p.x_) - p.y_ - 1.32};
}

Point next_point_Newton(Point p)
{
    Point f = f_Newton(p);
    double detJ = 1.0 + sin(p.y_) * cos(p.x_);
    double x_new = p.x_ - (-1 * f.x_ + sin(p.y_) * f.y_) / detJ;
    double y_new = p.y_ - (-cos(p.x_) * f.x_ - f.y_) / detJ;
    return Point{x_new, y_new};
}

double err(const Point& p)
{
    return std::max(fabs(cos(p.y_) - p.x_ + 0.85), fabs(sin(p.x_) - p.y_ - 1.32));
}

std::ostream& operator<<(std::ostream& out, const Point& p) {
    out << "Point{" << p.x_ << ", " << p.y_ << "}  err = " << err(p) << "\n";
    return out;
}

int main()
{
    std::ofstream fout;
    fout.open("out.txt");

    fout << "------MPI------\n";

    Point p{1.8, -0.4};
    fout << p;

    while(err(p) > ERR) {
        p = next_point_MPI(std::move(p));
        fout << p;
    }

    fout << "-----Newton-----\n";

    p = Point{1.8, -0.4};
    fout << p;

    while(err(p) > ERR) {
        p = next_point_Newton(std::move(p));
        fout << p;
    }

    fout.close();
    return 0;
}
