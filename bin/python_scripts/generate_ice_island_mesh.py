import struct
import math

nums = [60, 60, 60]
lengths = [20.0, 20.0, 20.0]

ground_z = 2.0
water_z = 15.0
ice_island_z = 17.0
ice_island_radius = 5.0
ground_temperature = 250.0
water_temperature = 278.0
ice_temperature = 270.0
air_temperature = 270.0



def generate_ice_island(file_name, nums, lengths, ground_num_z, water_num_z, ice_island_num_z, ice_island_radius_num,
	temperatures_file_name = None, temperature_by_param = None):
	dims = len(nums)
	steps = []
	for i in range(len(nums)):
		steps.append(lengths[i] / nums[i])
	num_of_nodes = reduce(lambda x, y: x*y, nums)
	file = open(file_name, 'wb')
	file.write(struct.pack("<%ui" % len(nums), *nums))
	file.write(struct.pack("<%ud" % len(lengths), *lengths))

	if (temperatures_file_name):
		temperatures_file = open(temperatures_file_name, 'wb')

	params = []
	temperatures = []
	center_x = lengths[0] / 2.0
	center_y = lengths[1] / 2.0
	for k in range(nums[2]):
		for j in range(nums[1]):
			for i in range(nums[0]):
				x = i * lengths[0] / nums[0] - center_x
				y = 0.0
				z = 0.0
				r = 0.0
				if dims == 3:
					y = j * lengths[1] / nums[1] - center_y
					z = k * lengths[2] / nums[2]
					r = math.sqrt(x*x + y*y)
				elif dims == 2:
					z = j * lengths[1] / nums[1]
					r = math.abs(x)

				if z < ground_z:
					params.append(0)
				elif z < water_z or (z < ice_island_z and r <= ice_island_radius):
					params.append(1)
				else:
					params.append(2)
				if (temperatures_file_name):
					if z < ground_z:
						temperatures.append(temperature_by_param[0])
					elif z < water_z and r > ice_island_radius:
						temperatures.append(temperature_by_param[1])
					elif z < ice_island_z and r <= ice_island_radius:
						temperatures.append(temperature_by_param[2])
					else:
						temperatures.append(temperature_by_param[3])


	file.write(struct.pack("<%ui" % len(params), *params))
	file.close()
	temperatures_file.write(struct.pack("<%ud" % len(temperatures), *temperatures))
	temperatures_file.close()



generate_ice_island("../meshes/ice_island3d.mesh", nums, lengths, ground_z, water_z, ice_island_z, ice_island_radius,
	"../meshes/ice_island3d_temperature.data", [ground_temperature, water_temperature, ice_temperature, air_temperature])
