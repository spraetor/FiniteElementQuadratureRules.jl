# Quadrature Rule Overview

This page lists all quadrature rules found in `rules/compact`.

Quality markers:
- `P`: positive weights
- `N`: some negative weights
- `I`: all quadrature points are strictly inside the reference element
- `B`: some quadrature points lie on the boundary
- `O`: some quadrature points lie outside the reference element

## Hexahedron

| Degree | Points | Quality | Orbits | Reference |
| ---: | ---: | :---: | :--- | :--- |
| 1 | 1 | PI | `[1, 0, 0, 0, 0, 0, 0]` | [`WV15`](#WV15) |
| 3 | 6 | PB | `[0, 1, 0, 0, 0, 0, 0]` | [`WV15`](#WV15) |
| 5 | 14 | PI | `[0, 1, 1, 0, 0, 0, 0]` | [`WV15`](#WV15) |
| 7 | 34 | PI | `[0, 1, 2, 1, 0, 0, 0]` | [`WV15`](#WV15) |
| 9 | 58 | PI | `[0, 1, 2, 1, 0, 1, 0]` | [`WV15`](#WV15) |
| 11 | 90 | PI | `[0, 1, 3, 1, 0, 2, 0]` | [`WV15`](#WV15) |

## Prism

| Degree | Points | Quality | Orbits | Reference |
| ---: | ---: | :---: | :--- | :--- |
| 1 | 1 | PI | `[1, 0, 0, 0, 0, 0]` | [`WV15`](#WV15) |
| 2 | 5 | PI | `[0, 1, 1, 0, 0, 0]` | [`WV15`](#WV15) |
| 3 | 8 | PI | `[0, 1, 0, 0, 1, 0]` | [`WV15`](#WV15) |
| 4 | 11 | PI | `[0, 1, 1, 1, 0, 0]` | [`WV15`](#WV15) |
| 5 | 16 | PI | `[1, 0, 1, 2, 0, 0]` | [`WV15`](#WV15) |
| 6 | 28 | PI | `[0, 2, 2, 1, 0, 1]` | [`WV15`](#WV15) |
| 7 | 35 | PI | `[0, 1, 1, 2, 1, 1]` | [`WV15`](#WV15) |
| 8 | 46 | PI | `[0, 2, 2, 4, 0, 1]` | [`WV15`](#WV15) |
| 9 | 60 | PI | `[1, 1, 1, 6, 1, 1]` | [`WV15`](#WV15) |
| 10 | 85 | PB | `[0, 2, 3, 5, 1, 3]` | [`WV15`](#WV15) |

## Pyramid

| Degree | Points | Quality | Orbits | Reference |
| ---: | ---: | :---: | :--- | :--- |
| 1 | 1 | PI | `[1, 0, 0, 0]` | [`JS21`](#JS21) |
| 1 | 1 | PI | `[1, 0, 0, 0]` | [`WV15`](#WV15) |
| 2 | 5 | PI | `[1, 0, 1, 0]` | [`JS21`](#JS21) |
| 2 | 5 | PI | `[1, 1, 0, 0]` | [`WV15`](#WV15) |
| 3 | 6 | PI | `[2, 0, 1, 0]` | [`JS21`](#JS21) |
| 3 | 6 | PI | `[2, 0, 1, 0]` | [`WV15`](#WV15) |
| 4 | 10 | PI | `[2, 1, 1, 0]` | [`JS21`](#JS21) |
| 4 | 10 | PI | `[2, 1, 1, 0]` | [`WV15`](#WV15) |
| 5 | 15 | PI | `[3, 1, 2, 0]` | [`JS21`](#JS21) |
| 5 | 15 | PI | `[3, 1, 2, 0]` | [`WV15`](#WV15) |
| 6 | 23 | PI | `[3, 2, 3, 0]` | [`JS21`](#JS21) |
| 6 | 24 | PB | `[4, 2, 3, 0]` | [`WV15`](#WV15) |
| 7 | 31 | PI | `[3, 2, 5, 0]` | [`JS21`](#JS21) |
| 7 | 31 | PI | `[3, 2, 5, 0]` | [`WV15`](#WV15) |
| 8 | 47 | PI | `[3, 3, 6, 1]` | [`WV15`](#WV15) |
| 8 | 47 | PI | `[3, 4, 5, 1]` | [`JS21`](#JS21) |
| 9 | 62 | PI | `[2, 4, 9, 1]` | [`WV15`](#WV15) |
| 9 | 62 | PI | `[2, 6, 5, 2]` | [`JS21`](#JS21) |
| 10 | 80 | PI | `[4, 6, 7, 3]` | [`JS21`](#JS21) |
| 10 | 83 | PI | `[3, 3, 9, 4]` | [`WV15`](#WV15) |
| 11 | 103 | PI | `[3, 8, 11, 3]` | [`JS21`](#JS21) |
| 12 | 127 | PI | `[3, 9, 12, 5]` | [`JS21`](#JS21) |
| 13 | 152 | PI | `[4, 8, 15, 7]` | [`JS21`](#JS21) |
| 14 | 184 | PI | `[4, 12, 17, 8]` | [`JS21`](#JS21) |
| 15 | 234 | PI | `[2, 13, 23, 11]` | [`JS21`](#JS21) |
| 16 | 285 | PI | `[1, 13, 22, 18]` | [`JS21`](#JS21) |
| 17 | 319 | PI | `[3, 13, 28, 19]` | [`JS21`](#JS21) |
| 18 | 357 | PI | `[1, 14, 31, 22]` | [`JS21`](#JS21) |
| 19 | 418 | PI | `[2, 17, 29, 29]` | [`JS21`](#JS21) |
| 20 | 489 | PI | `[1, 19, 31, 36]` | [`JS21`](#JS21) |

## Quadrilateral

| Degree | Points | Quality | Orbits | Reference |
| ---: | ---: | :---: | :--- | :--- |
| 1 | 1 | PI | `[1, 0, 0, 0]` | [`WV15`](#WV15) |
| 3 | 4 | PI | `[0, 0, 1, 0]` | [`WV15`](#WV15) |
| 5 | 8 | PI | `[0, 1, 1, 0]` | [`WV15`](#WV15) |
| 7 | 12 | PI | `[0, 1, 2, 0]` | [`WV15`](#WV15) |
| 9 | 20 | PI | `[0, 1, 2, 1]` | [`WV15`](#WV15) |
| 11 | 28 | PI | `[0, 1, 2, 2]` | [`WV15`](#WV15) |
| 13 | 37 | PI | `[1, 2, 3, 2]` | [`WV15`](#WV15) |
| 15 | 48 | PI | `[0, 2, 2, 4]` | [`WV15`](#WV15) |
| 17 | 60 | PI | `[0, 2, 5, 4]` | [`WV15`](#WV15) |
| 19 | 72 | PI | `[0, 3, 3, 6]` | [`WV15`](#WV15) |
| 21 | 85 | PI | `[1, 2, 5, 7]` | [`WV15`](#WV15) |

## Tetrahedron

| Degree | Points | Quality | Orbits | Reference |
| ---: | ---: | :---: | :--- | :--- |
| 1 | 1 | PI | `[1, 0, 0, 0, 0]` | [`JS21`](#JS21) |
| 1 | 1 | PI | `[1, 0, 0, 0, 0]` | [`WV15`](#WV15) |
| 2 | 4 | PI | `[0, 1, 0, 0, 0]` | [`JS21`](#JS21) |
| 2 | 4 | PI | `[0, 1, 0, 0, 0]` | [`WV15`](#WV15) |
| 3 | 8 | PI | `[0, 2, 0, 0, 0]` | [`JS21`](#JS21) |
| 3 | 8 | PI | `[0, 2, 0, 0, 0]` | [`WV15`](#WV15) |
| 4 | 14 | PI | `[0, 2, 1, 0, 0]` | [`CCGV22`](#CCGV22) |
| 4 | 14 | PI | `[0, 2, 1, 0, 0]` | [`JS21`](#JS21) |
| 5 | 14 | PI | `[0, 2, 1, 0, 0]` | [`CCGV22`](#CCGV22) |
| 5 | 14 | PI | `[0, 2, 1, 0, 0]` | [`JS21`](#JS21) |
| 5 | 14 | PI | `[0, 2, 1, 0, 0]` | [`WV15`](#WV15) |
| 6 | 24 | PI | `[0, 3, 0, 1, 0]` | [`CCGV22`](#CCGV22) |
| 6 | 24 | PI | `[0, 3, 0, 1, 0]` | [`JS21`](#JS21) |
| 6 | 24 | PI | `[0, 3, 0, 1, 0]` | [`WV15`](#WV15) |
| 7 | 35 | PI | `[1, 1, 1, 2, 0]` | [`CCGV22`](#CCGV22) |
| 7 | 35 | PI | `[1, 1, 1, 2, 0]` | [`JS21`](#JS21) |
| 7 | 35 | PI | `[1, 1, 1, 2, 0]` | [`WV15`](#WV15) |
| 8 | 46 | PI | `[0, 4, 1, 2, 0]` | [`CCGV22`](#CCGV22) |
| 8 | 46 | PI | `[0, 4, 1, 2, 0]` | [`JS21`](#JS21) |
| 8 | 46 | PI | `[0, 4, 1, 2, 0]` | [`WV15`](#WV15) |
| 8 | 46 | PI | `[0, 4, 1, 2, 0]` | [`ZCL09`](#ZCL09) |
| 9 | 59 | PI | `[1, 4, 1, 3, 0]` | [`CCGV22`](#CCGV22) |
| 9 | 59 | PI | `[1, 4, 1, 3, 0]` | [`JS21`](#JS21) |
| 9 | 59 | PI | `[1, 4, 1, 3, 0]` | [`WV15`](#WV15) |
| 10 | 79 | PI | `[1, 3, 1, 5, 0]` | [`CCGV22`](#CCGV22) |
| 10 | 81 | PI | `[1, 2, 0, 6, 0]` | [`WV15`](#WV15) |
| 10 | 81 | PI | `[1, 2, 2, 5, 0]` | [`JS21`](#JS21) |
| 11 | 98 | PI | `[0, 5, 1, 4, 1]` | [`CCGV22`](#CCGV22) |
| 11 | 110 | PI | `[0, 2, 3, 5, 1]` | [`JS21`](#JS21) |
| 12 | 123 | PI | `[1, 5, 1, 6, 1]` | [`CCGV22`](#CCGV22) |
| 12 | 168 | PI | `[0, 3, 2, 4, 4]` | [`JS21`](#JS21) |
| 13 | 145 | PI | `[1, 3, 2, 8, 1]` | [`CCGV22`](#CCGV22) |
| 13 | 172 | PI | `[0, 4, 2, 6, 3]` | [`JS21`](#JS21) |
| 14 | 175 | PI | `[1, 6, 1, 10, 1]` | [`CCGV22`](#CCGV22) |
| 14 | 204 | PI | `[0, 6, 4, 5, 4]` | [`JS21`](#JS21) |
| 14 | 236 | PI | `[0, 5, 0, 4, 7]` | [`ZCL09`](#ZCL09) |
| 15 | 209 | PI | `[1, 4, 2, 11, 2]` | [`CCGV22`](#CCGV22) |
| 15 | 264 | PI | `[0, 3, 2, 6, 7]` | [`JS21`](#JS21) |
| 16 | 248 | PI | `[0, 8, 2, 11, 3]` | [`CCGV22`](#CCGV22) |
| 16 | 304 | PI | `[0, 4, 2, 7, 8]` | [`JS21`](#JS21) |
| 17 | 284 | PI | `[0, 8, 2, 14, 3]` | [`CCGV22`](#CCGV22) |
| 17 | 364 | PI | `[0, 4, 4, 9, 9]` | [`JS21`](#JS21) |
| 18 | 343 | PI | `[1, 6, 1, 18, 4]` | [`CCGV22`](#CCGV22) |
| 18 | 436 | PI | `[0, 7, 8, 10, 10]` | [`JS21`](#JS21) |
| 19 | 383 | PI | `[1, 7, 3, 18, 5]` | [`CCGV22`](#CCGV22) |
| 19 | 487 | PI | `[1, 3, 1, 13, 13]` | [`JS21`](#JS21) |
| 20 | 441 | PI | `[1, 8, 4, 20, 6]` | [`CCGV22`](#CCGV22) |
| 20 | 552 | PI | `[0, 6, 2, 13, 15]` | [`JS21`](#JS21) |

## Triangle

| Degree | Points | Quality | Orbits | Reference |
| ---: | ---: | :---: | :--- | :--- |
| 1 | 1 | PI | `[1, 0, 0, 0]` | [`Pap15`](#Pap15) |
| 1 | 1 | PI | `[1, 0, 0]` | [`WV15`](#WV15) |
| 2 | 3 | PB | `[0, 1, 0, 0]` | [`LJ75`](#LJ75) |
| 2 | 3 | PB | `[0, 1, 0, 0]` | [`Pap15`](#Pap15) |
| 2 | 3 | PI | `[0, 1, 0]` | [`WV15`](#WV15) |
| 2 | 4 | NI | `[1, 1, 0, 0]` | [`LJ75`](#LJ75) |
| 2 | 4 | PB | `[1, 1, 0, 0]` | [`LJ75`](#LJ75) |
| 2 | 7 | PB | `[1, 2, 0, 0]` | [`Pap15`](#Pap15) |
| 2 | 7 | PI | `[1, 2, 0, 0]` | [`Pap15`](#Pap15) |
| 2 | 9 | PB | `[0, 3, 0, 0]` | [`Pap15`](#Pap15) |
| 3 | 4 | NI | `[1, 1, 0, 0]` | [`Pap15`](#Pap15) |
| 3 | 6 | PI | `[0, 0, 1, 0]` | [`Pap15`](#Pap15) |
| 4 | 6 | PI | `[0, 2, 0, 0]` | [`CCGV22`](#CCGV22) |
| 4 | 6 | PI | `[0, 2, 0, 0]` | [`LJ75`](#LJ75) |
| 4 | 6 | PI | `[0, 2, 0, 0]` | [`Pap15`](#Pap15) |
| 4 | 6 | PI | `[0, 2, 0]` | [`WV15`](#WV15) |
| 4 | 10 | NB | `[1, 1, 1, 0]` | [`LJ75`](#LJ75) |
| 5 | 7 | PI | `[1, 2, 0, 0]` | [`Pap15`](#Pap15) |
| 5 | 7 | PI | `[1, 2, 0]` | [`CCGV22`](#CCGV22) |
| 5 | 7 | PI | `[1, 2, 0]` | [`WV15`](#WV15) |
| 5 | 10 | PB | `[1, 3, 0, 0]` | [`Pap15`](#Pap15) |
| 6 | 12 | PI | `[0, 2, 1, 0]` | [`Pap15`](#Pap15) |
| 6 | 12 | PI | `[0, 2, 1, 0]` | [`Pap15`](#Pap15) |
| 6 | 12 | PI | `[0, 2, 1]` | [`CCGV22`](#CCGV22) |
| 6 | 12 | PI | `[0, 2, 1]` | [`WV15`](#WV15) |
| 6 | 16 | NB | `[1, 3, 1, 0]` | [`Pap15`](#Pap15) |
| 7 | 12 | PI | `[0, 0, 0, 0, 4]` | [`Gat88`](#Gat88) |
| 7 | 13 | NI | `[1, 2, 1, 0]` | [`Pap15`](#Pap15) |
| 7 | 15 | PI | `[0, 1, 2, 0]` | [`Pap15`](#Pap15) |
| 7 | 15 | PI | `[0, 1, 2]` | [`CCGV22`](#CCGV22) |
| 7 | 15 | PI | `[0, 3, 1]` | [`WV15`](#WV15) |
| 8 | 16 | NI | `[1, 3, 1, 0]` | [`Pap15`](#Pap15) |
| 8 | 16 | PI | `[1, 3, 1]` | [`CCGV22`](#CCGV22) |
| 8 | 16 | PI | `[1, 3, 1]` | [`WV15`](#WV15) |
| 8 | 16 | PI | `[1, 3, 1]` | [`ZCL09`](#ZCL09) |
| 9 | 19 | PI | `[1, 4, 1, 0]` | [`Pap15`](#Pap15) |
| 9 | 19 | PI | `[1, 4, 1]` | [`CCGV22`](#CCGV22) |
| 9 | 19 | PI | `[1, 4, 1]` | [`WV15`](#WV15) |
| 10 | 24 | PI | `[0, 0, 0, 0, 8]` | [`Pap16`](#Pap16) |
| 10 | 24 | PI | `[0, 0, 0, 10, 0, 4]` | [`XG10`](#XG10) |
| 10 | 24 | PO | `[0, 4, 2, 0]` | [`Pap15`](#Pap15) |
| 10 | 25 | PI | `[1, 2, 3, 0]` | [`Pap15`](#Pap15) |
| 10 | 25 | PI | `[1, 2, 3]` | [`CCGV22`](#CCGV22) |
| 10 | 25 | PI | `[1, 2, 3]` | [`WV15`](#WV15) |
| 11 | 27 | PI | `[0, 0, 0, 0, 9]` | [`Pap16`](#Pap16) |
| 11 | 27 | PO | `[0, 5, 2, 0]` | [`Pap15`](#Pap15) |
| 11 | 28 | NI | `[1, 3, 3, 0]` | [`Pap15`](#Pap15) |
| 11 | 28 | PI | `[1, 5, 2]` | [`CCGV22`](#CCGV22) |
| 11 | 28 | PI | `[1, 5, 2]` | [`WV15`](#WV15) |
| 11 | 30 | PI | `[0, 2, 4, 0]` | [`Pap15`](#Pap15) |
| 12 | 33 | NI | `[0, 5, 3, 0]` | [`Pap15`](#Pap15) |
| 12 | 33 | PI | `[0, 5, 3]` | [`CCGV22`](#CCGV22) |
| 12 | 33 | PI | `[0, 5, 3]` | [`WV15`](#WV15) |
| 13 | 36 | PI | `[0, 0, 0, 0, 12]` | [`Pap16`](#Pap16) |
| 13 | 37 | PI | `[1, 4, 4, 0]` | [`Pap15`](#Pap15) |
| 13 | 37 | PI | `[1, 4, 4]` | [`CCGV22`](#CCGV22) |
| 13 | 37 | PI | `[1, 4, 4]` | [`WV15`](#WV15) |
| 14 | 40 | PO | `[1, 0, 0, 0, 13]` | [`Pap16`](#Pap16) |
| 14 | 42 | PI | `[0, 6, 4, 0]` | [`Pap15`](#Pap15) |
| 14 | 42 | PI | `[0, 6, 4]` | [`CCGV22`](#CCGV22) |
| 14 | 42 | PI | `[0, 6, 4]` | [`WV15`](#WV15) |
| 14 | 46 | PI | `[1, 7, 4]` | [`ZCL09`](#ZCL09) |
| 15 | 46 | PI | `[1, 0, 0, 0, 15]` | [`Pap16`](#Pap16) |
| 15 | 48 | NI | `[0, 6, 5, 0]` | [`Pap16`](#Pap16) |
| 15 | 49 | PI | `[1, 4, 6]` | [`CCGV22`](#CCGV22) |
| 15 | 49 | PI | `[1, 6, 5]` | [`WV15`](#WV15) |
| 16 | 51 | PO | `[0, 0, 0, 0, 17]` | [`Pap16`](#Pap16) |
| 16 | 52 | PI | `[1, 0, 0, 0, 17]` | [`Pap16`](#Pap16) |
| 16 | 55 | PB | `[1, 6, 6]` | [`WV15`](#WV15) |
| 16 | 55 | PI | `[1, 4, 7]` | [`CCGV22`](#CCGV22) |
| 17 | 57 | PI | `[0, 0, 0, 0, 19]` | [`Pap16`](#Pap16) |
| 17 | 58 | PO | `[1, 7, 6, 0]` | [`Pap16`](#Pap16) |
| 17 | 60 | PI | `[0, 6, 7]` | [`CCGV22`](#CCGV22) |
| 17 | 60 | PI | `[0, 6, 7]` | [`WV15`](#WV15) |
| 18 | 64 | NI | `[1, 0, 0, 0, 21]` | [`Pap16`](#Pap16) |
| 18 | 66 | PI | `[0, 0, 0, 0, 22]` | [`Pap16`](#Pap16) |
| 18 | 67 | PI | `[1, 6, 8]` | [`CCGV22`](#CCGV22) |
| 18 | 67 | PI | `[1, 6, 8]` | [`WV15`](#WV15) |
| 19 | 70 | PI | `[1, 0, 0, 0, 23]` | [`Pap16`](#Pap16) |
| 19 | 73 | PI | `[1, 6, 9]` | [`CCGV22`](#CCGV22) |
| 19 | 73 | PI | `[1, 8, 8]` | [`WV15`](#WV15) |
| 20 | 78 | PI | `[0, 0, 0, 0, 26]` | [`Pap16`](#Pap16) |
| 20 | 79 | PI | `[1, 8, 9]` | [`CCGV22`](#CCGV22) |
| 20 | 79 | PI | `[1, 8, 9]` | [`WV15`](#WV15) |
| 20 | 88 | PI | `[1, 5, 12]` | [`ZCL09`](#ZCL09) |
| 21 | 85 | PI | `[1, 0, 0, 0, 28]` | [`Pap16`](#Pap16) |
| 21 | 87 | PI | `[0, 9, 10, 0]` | [`Pap16`](#Pap16) |
| 22 | 93 | PI | `[0, 0, 0, 0, 31]` | [`Pap16`](#Pap16) |
| 22 | 96 | PI | `[0, 10, 11, 0]` | [`Pap16`](#Pap16) |
| 23 | 100 | PI | `[1, 0, 0, 0, 33]` | [`Pap16`](#Pap16) |
| 23 | 102 | PI | `[0, 10, 12, 0]` | [`Pap16`](#Pap16) |
| 24 | 109 | PI | `[1, 0, 0, 0, 36]` | [`Pap16`](#Pap16) |
| 24 | 112 | PI | `[1, 9, 14, 0]` | [`Pap16`](#Pap16) |
| 25 | 117 | PI | `[0, 0, 0, 0, 39]` | [`Pap16`](#Pap16) |
| 25 | 118 | PO | `[1, 11, 14, 0]` | [`Pap16`](#Pap16) |
| 25 | 120 | NI | `[0, 12, 14, 0]` | [`Pap16`](#Pap16) |
| 30 | 168 | PI | `[0, 0, 0, 0, 56]` | [`XG10`](#XG10) |
| 50 | 453 | PI | `[0, 17, 67]` | [`XG10`](#XG10) |

## References

### CCGV22

G. Chuluunbaatar, O. Chuluunbaatar, A. A. Gusev, and S. I. Vinitsky.
Pi-type fully symmetric quadrature rules on the 3-,$\ldots$, 6-simplexes.
*Computer & Mathematics with Applications*, 124:89–97, 2022.

### Gat88

K. Gatermann.
The construction of symmetric cubature formulas for the square and the triangle.
*Computing*, 40:229–240, 1988.

### JS21

J. Jaśkowiec and N. Sukumar.
High-order symmetric cubature rules for tetrahedra and pyramids.
*International Journal for Numerical Methods in Engineering*, 122(1):148–171, 2021.

### LJ75

J. N. Lyness and D. Jespersen.
Moderate degree symmetric quadrature rules for the triangle.
*J. Inst. Math. Appl.*, 15(1):19–32, 1975.

### Pap15

S.-A. Papanicolopulos.
Computation of moderate-degree fully-symmetric cubature rules on the triangle using symmetric polynomials and algebraic solving.
*Computers & Mathematics with Applications*, 69:650–666, 2015.

### Pap16

S.-A. Papanicolopulos.
New fully symmetric and rotationally symmetric cubature rules on the triangle using minimal orthonormal bases.
*Journal of Computational and Applied Mathematics*, 294:39–48, 2016.

### WV15

F. D. Witherden and P. E. Vincent.
On the identification of symmetric quadrature rules for finite element methods.
*Computers & Mathematics with Applications*, 69, 2015.

### XG10

H. Xiao and Z. Gimbutas.
A numerical algorithm for the construction of efficient quadrature rules in two and higher dimensions.
*Computers & Mathematics with Applications*, 59(2):663–676, 2010.

### ZCL09

L. Zhang, T. Cui, and H. Liu.
A set of symmetric quadrature rules on triangles and tetrahedra.
*Journal of Computational Mathematics*, 27(1):89–96, 2009.
