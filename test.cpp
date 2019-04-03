#include <iostream>
#include <complex>

class plan_s {
    public:
        virtual void f() = 0;
        virtual ~plan_s() { std::cout << "deleting parent\n"; }
        virtual void execute() = 0;
};

class plan_s1 : public plan_s {
    public:
        virtual void f() { std::cout << "1\n"; }
        virtual ~plan_s1() { std::cout << "deleting 1\n"; }
        virtual void execute() { std::cout << "executing plan 1" << std::endl; }
};

class plan_s2 : public plan_s {
    public:
        virtual void f() { std::cout << "2\n"; }
        virtual ~plan_s2() { std::cout << "deleting 2\n"; }
        virtual void execute() { std::cout << "executing plan 2" << std::endl; }
};

typedef plan_s* plan;

plan make_plan1() {
    plan p = new plan_s1;
    return p;
}

plan make_plan2() {
    plan p = new plan_s2;
    return p;
}

void destroy_plan(plan p) {
    if (p) delete p;
}

void execute(plan p) {
    p->execute();
}

template<int N> class node;

template<>
class node<1> {
    int x;
};

template<int N>
class node {
    node<1> ns[N];
};

void cx(std::complex<double> x) {
    std::cout << x.real() << std::endl;
}

int main() {
    plan p1 = make_plan1();
    p1->f();
    plan p2 = make_plan2();
    p2->f();

    execute(p1);
    execute(p2);

    destroy_plan(p1);
    destroy_plan(p2);

    node<2> n;

    double x = 2;
    cx(x);
}
