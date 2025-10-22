Ni clang:

```shell
iota_y_zip.cpp:10:33: error: no member named 'views' in namespace 'std'
10         | for (auto [i, valor] : std::views::zip(std::views::iota(0), huevo)) {
           | ~~~~~^
        iota_y_zip.cpp:10:49: error: no member named 'views' in namespace 'std'
        10 | for (auto [i, valor] : std::views::zip(std::views::iota(0), huevo)) {
           | ~~~~~^
                2 errors generated.
```

Ni g++:

```shell
iota_y_zip.cpp: In function ‘int main()’:
iota_y_zip.cpp:10:33: error: ‘std::views’ has not been declared
   10 |     for (auto [i, valor] : std::views::zip(std::views::iota(0), huevo)) {
      |                                 ^~~~~
iota_y_zip.cpp:10:49: error: ‘std::views’ has not been declared
   10 |     for (auto [i, valor] : std::views::zip(std::views::iota(0), huevo)) {
      |                                                 ^~~~~
```

detectan a `std::views`, las versiones de los compiladores son las siguientes:

gcc: `g++ (Ubuntu 13.3.0-6ubuntu2~24.04) 13.3.0`
clang: `Ubuntu clang version 18.1.3 (1ubuntu1)`
