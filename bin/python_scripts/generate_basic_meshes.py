import struct
import math
from functools import reduce


def generate_zeroes(file_name, nums, lengths):
	dims = len(nums)
	num_of_nodes = reduce(lambda x, y: x*y, nums)
	file = open(file_name, 'wb')
	file.write(struct.pack("<%ui" % len(nums), *nums))
	file.write(struct.pack("<%ud" % len(lengths), *lengths))
	params = [0] * num_of_nodes
	file.write(struct.pack("<%ui" % len(params), *params))
	file.close()

def generate_half_zeroes_half_ones(file_name, nums, lengths):
	dims = len(nums)
	num_of_nodes = reduce(lambda x, y: x*y, nums)
	file = open(file_name, 'wb')
	file.write(struct.pack("<%ui" % len(nums), *nums))
	file.write(struct.pack("<%ud" % len(lengths), *lengths))
	params = [0] * (num_of_nodes//2) + [1] * (num_of_nodes//2)
	file.write(struct.pack("<%ui" % len(params), *params))
	file.close()


def generate_box_2d(file_name, nums, lengths, size):
	dims = len(nums)
	num_of_nodes = reduce(lambda x, y: x*y, nums)
	file = open(file_name, 'wb')
	file.write(struct.pack("<%ui" % len(nums), *nums))
	file.write(struct.pack("<%ud" % len(lengths), *lengths))

	params = []
	center_x = lengths[0] / 2.0
	center_y = lengths[1] / 2.0
	for i in range(nums[0]):
		for j in range(nums[1]):
			x = i * lengths[0] / nums[0] - center_x
			y = j * lengths[1] / nums[1] - center_y
			if math.fabs(x) < size/2 and math.fabs(y) < size/2:
				params.append(1)
			else:
				params.append(0)

	file.write(struct.pack("<%ui" % len(params), *params))
	file.close()


def generate_box_3d(file_name, nums, lengths, size):
	dims = len(nums)
	num_of_nodes = reduce(lambda x, y: x*y, nums)
	file = open(file_name, 'wb')
	file.write(struct.pack("<%ui" % len(nums), *nums))
	file.write(struct.pack("<%ud" % len(lengths), *lengths))

	params = []
	center_x = lengths[0] / 2.0
	center_y = lengths[1] / 2.0
	center_z = lengths[2] / 2.0
	for i in range(nums[0]):
		for j in range(nums[1]):
			for k in range(nums[2]):
				x = i * lengths[0] / nums[0] - center_x
				y = j * lengths[1] / nums[1] - center_y
				z = k * lengths[2] / nums[2] - center_z
				if math.fabs(x) < size/2 and math.fabs(y) < size/2 and math.fabs(z) < size/2:
					params.append(1)
				else:
					params.append(0)

	file.write(struct.pack("<%ui" % len(params), *params))
	file.close()



def generate_ball_2d(file_name, nums, lengths, size):
	file = open(file_name, 'wb')
	file.write(struct.pack("<%ui" % len(nums), *nums))
	file.write(struct.pack("<%ud" % len(lengths), *lengths))

	params = []
	center_x = lengths[0] / 2.0
	center_y = lengths[1] / 2.0
	for i in range(nums[0]):
		for j in range(nums[1]):
			x = i * lengths[0] / nums[0] - center_x
			y = j * lengths[1] / nums[1] - center_y
			if x*x + y*y < size*size/4:
				params.append(1)
			else:
				params.append(0)

	file.write(struct.pack("<%ui" % len(params), *params))
	file.close()


def generate_ball_2d_ascii(file_name, nums, lengths, size):
	file = open(file_name, 'w')
	file.write(str(nums[0]) + " " + str(nums[1]) + "\n")
	file.write(str(lengths[0]) + " " + str(lengths[1]) + "\n")

	center_x = lengths[0] / 2.0
	center_y = lengths[1] / 2.0
	for i in range(nums[0]):
		for j in range(nums[1]):
			x = i * lengths[0] / nums[0] - center_x
			y = j * lengths[1] / nums[1] - center_y
			if x*x + y*y < size*size/4:
				file.write("1 ")
			else:
				file.write("0 ")
		file.write("\n")
	file.close()


def generate_ball_field_2d(file_name, nums, lengths, size):
	file = open(file_name, 'w')
	for i in range(nums[0]):
		for j in range(nums[1]):
			x = i * lengths[0] / nums[0]
			y = j * lengths[1] / nums[1]
			if y < 3:
				file.write(str(i) + " " + str(j) + "\n")
	file.close()


def generate_ball_3d(file_name, nums, lengths, size):
	file = open(file_name, 'wb')
	file.write(struct.pack("<%ui" % len(nums), *nums))
	file.write(struct.pack("<%ud" % len(lengths), *lengths))

	params = []
	center_x = lengths[0] / 2.0
	center_y = lengths[1] / 2.0
	center_z = lengths[2] / 2.0
	for i in range(nums[0]):
		for j in range(nums[1]):
			for k in range(nums[2]):
				x = i * lengths[0] / nums[0] - center_x
				y = j * lengths[1] / nums[1] - center_y
				z = k * lengths[2] / nums[2] - center_z
				if x*x + y*y + z*z < size*size/4:
					params.append(1)
				else:
					params.append(0)

	file.write(struct.pack("<%ui" % len(params), *params))
	file.close()



def generate_linear_temperatures(file_name, nums, axis, start_temperature, finish_temperature):
	dims = len(nums)
	file = open(file_name, 'wb')
	#file.write(struct.pack("<%ui" % len(nums), *nums))
	params = []

	if dims == 2:
		for i in range(nums[0]):
			for j in range(nums[1]):
				temperature_step = (finish_temperature - start_temperature) / nums[axis]
				index = (i, j)[axis]
				params.append(start_temperature + temperature_step * index)
	if dims == 3:
		for i in range(nums[0]):
			for j in range(nums[1]):
				for k in range(nums[2]):
					temperature_step = (finish_temperature - start_temperature) / float(nums[axis])
					index = (i, j, k)[axis]
					params.append(start_temperature + temperature_step * index)

	file.write(struct.pack("<%uf" % len(params), *params))
	file.close()

mode = "TEST"

if mode == "TEST":
	generate_ball_2d_ascii("../meshes/ball2d_512.tmesh", (512, 512), (10.0, 10.0), 4.0)
	generate_ball_field_2d("../meshes/ball2d_512.tfield", (512, 512), (10.0, 10.0), 4.0)
	

	#generate_ball_2d("../meshes/ball2d_1024.bmesh", (1024, 1024), (10.0, 10.0), 4.0)
	#generate_ball_3d("../meshes/ball3d_128.bmesh", (128, 128, 128), (10.0, 10.0, 10.0), 4.0)
	#generate_ball_3d("../meshes/ball3d_256.bmesh", (256, 256, 256), (10.0, 10.0, 10.0), 4.0)


