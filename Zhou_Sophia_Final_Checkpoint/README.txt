Link: https://youtu.be/UGgISMncI8I 

What effects you implemented, and where they are in the video. (1 point)
- Soft shadows and texture mapping. You can see them on the floor, the rock/"chair", and in the sky(box).

What was easy to get working? (3 points)
- Soft shadows were pretty straightforward.

What was hard to get working, or didn’t work at all? (3 points)
- It took me a while to understand and iron out texture mapping bugs. Cylinder intersection was also difficult.

If you did it over again, what would you do different? (3 points)
- I would write cleaner code and structure the repository in a more effective way. 
- I might also have chosen a different action.
- And added a different scripted camera motion.

Instructions to build and run the code:
Add the Eigen directory inside of the previz directory. In the previz directory, type the following commands into the terminal:

make -B 
./previz 0 1

When the program says it has rendered one frame, check the frames directory and the rendered frame should be in there called frame.0000.ppm.