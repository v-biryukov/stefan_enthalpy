import struct
import math

nums = [100, 50]
lengths = [20.0, 10.0]

width = lengths[0]
height = lengths[1]
# Proportions:
ground_percent = 0.1
water_percent = 0.6
ice_height_percent = 0.75
air_percent = 1.0 - ground_percent - water_percent

ice_width_percent = 0.7

def generate_ice_island(file_name):
	file = open(file_name, 'w')
	file.write(str(nums[0]) + " " + str(nums[1]) + "\n")
	file.write(str(lengths[0]) + " " + str(lengths[1]) + "\n")

	params = []
	temperatures = []
	center_x = lengths[0] / 2.0
	center_y = lengths[1] / 2.0
	for j in range(nums[1]):
		for i in range(nums[0]):
			x = i * lengths[0] / nums[0] - center_x
			y = j * lengths[1] / nums[1]
			r = abs(x)

			if y < ground_percent * height:
				file.write("0 ")
			elif y < (ground_percent + water_percent) * height and r > width * ice_width_percent / 2:
				file.write("1 ")
			elif y < (ground_percent + water_percent) * height and r <= width * ice_width_percent / 2:
				file.write("2 ")
			elif y < (ground_percent + ice_height_percent) * height and r <= width * ice_width_percent / 2:
				file.write("3 ")
			else:
				file.write("4 ")
		file.write("\n")
	file.close()



generate_ice_island("./meshes/ice_island2d.mesh")
