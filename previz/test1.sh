#!/bin/bash

ffmpeg -r 30 -f image2 -s 640x480 -start_number 1 -i frame.%04d.ppm -vframes 1000 -vcodec libx264 -crf 25 -pix_fmt yuv420p movie2.mp4