# DEM-slope-pathfinder

This repository contains a python script I have written for a course at the Technical University of Dresden. 
This was my first time creating an algorithm from scratch. 
Although the solution still contains errors and can of course be improved and be written more elegant it was a good learning experience.

## The problem
The main idea was to model serpentines for keeping a certain profile below a defined gradient. 
First, a line would be drawn in a DEM as a proposed route. Along this line for each segment (spaced using DEM resolution) the gradient will be calculated.
If this gradient is larger than a set gradient the line is cut and a different route is calculated to reach the next point. 

## Proposed solution 
A filter is used which moves along the path. At each point of the DEM along the path the gradient to all neighbours in front are calculated. 
If the gradient is lower than the threshold this position is taken, the following constraints have to be taken into account:

1. the angle is lower than the threshold
2. the movement is towards the end-point
3. If more angles are lower take the lowest AND
4. take the one point towards the shortest distance to the end point.

Store all these points into a dictionary. 

## Difficulties
At first, the algorithm was developed and tested for a south to north facing case.
However, paths can also move from north to south and from west to east and vice versa. 
To calculate, for example, west to east moving route the numpy array is flipped to be orientated as south to north facing processed and flipped back and stored.

Moreover, starting at an outer edges caused an out of bound error therefore a rim with 'no data value' is added.

## Task list

- [x] Provide a readme and publish algorithm 
- [ ] Write down improvements and bugs in the section 'Improvements and bugs' 
- [ ] Make these changes

## Improvements and bugs

-- to be added --
