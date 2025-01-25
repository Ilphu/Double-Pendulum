# Double Pendulums and Heatmap

![image](heatmap.png)

## Description

A simple C++ prgoram designed to find and display stable double pendulum orbits.

## Compilation and Use:

To compile this program, run `cmake -B build` then `cmake --build build` in a shell of your choice. Run `./build/bin/dp -r` to render a stable pendulum orbit or run `./build/bin/dp -hm 100` to render a 100x100 heatmap where 100 can be changed to whatever you like.

## Credits:

All images created in this project are created using the `stbi` libraries.

Major inspo for this project comes form Joe Liang, highly recommend you go check out his stuff. This is ensentially the same thing as what he did however written in C++. All of the code not in the `stbi` files is mine, I intentionally didn't look at his repo. I also got insporoation from many others who have doen similar things to, MIT, and many youtube videos are what inspired me. 