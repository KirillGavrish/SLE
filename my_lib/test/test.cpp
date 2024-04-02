template <typename T>
void f1(T &a, T &b) {a = a + b;};

template <typename T, template <typename P> typename ftype>
void functor(T &a = 0};

int main()
{
    int a = 1;
    functor<f>(a);
}
