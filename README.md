# MCL
Monte Carlo Library

## Summary

Implementations of Monte Carlo integration schemes and samplers. To be used with [The Stochastic Simulation Algorithm Library (SSAL)](https://github.com/davidwarne/SSAL)

## Developer
David J. Warne (david.warne@qut.edu.au), School of Mathematical Sciences, Science and Engineering Faculty, Queensland University of Technology.

## License

MCL: Monte Carlo Library
Copyright (C) 2018  David J. Warne

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.


## Current Features

1. Monte Carlo integration routines for computing expectations of [SSAL](https://github.com/davidwarne/SSAL) based models, including multilevel Monte Carlo methods.
2. Various samplers for approximate Bayesian compuation (ABC), rejection sampling, Markov chain Monte Carlo, sequential Monte Carlo, and multilevel Monte Carlo 
3. Pseudo marginal methods

The library has been designed with performance and usability in mind. 

## Future development plans

Many of these rountines I use in my own research. Therefore, the current features and level of performance etc will generally be updated based on my own usage of the library. I am planning parallel (multi-threaded and distributed) version of many of the rountines.

If you find this library useful, or have any suggestions for improvement, I am always happy to be contacted.
