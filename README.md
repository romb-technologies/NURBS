# NURBS
A C++ library for 2D NURBS curves.
The project includes an example program to demonstrate the library's features.

## Features

- Define and evaluate NURBS curves of any order
- Manipulate control points, weights and knot vectors
- Evaluate derivatives, tangents, point projections, curve length, curvature etc.
- Control point and knot insertion, knot refinement, curve splitting
- Piecewise linear representation for drawing
- And more :smile:

The example program shows how the library can be integrated with Qt, as well as demonstrating 
most of the library's features interactively.

## Screenshots
![image showing the UI of the example program, with two curves](img/nurbs.png)

## Dependencies
- The library requires Eigen 3.3 or higher (https://eigen.tuxfamily.org/)
- The example program is made using Qt5 (https://www.qt.io/download)

## Installation
```
git clone https://github.com/romb-technologies/NURBS
mkdir NURBS/build
cd NURBS/build
cmake ..
make
make install
```

## Documentation
Doxygen is required to build the documentation locally.

1. In a terminal inside the project folder, run:
        ```
        doxygen
        ```
2. Navigate to `docs/html` and open `index.html`

## Licence and credits
Apache License Version 2.0

The example program uses [QCustomPlot](https://www.qcustomplot.com/)

---

Created at Romb Technologies (https://romb-technologies.hr/)
