import struct
import math

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
	params = [0] * (num_of_nodes/2) + [1] * (num_of_nodes/2)
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



generate_zeroes("../meshes/zeroes3d.mesh", (32, 32, 32), (10.0, 10.0, 10.0))
generate_zeroes("../meshes/zeroes2d.mesh", (64, 64), (10.0, 10.0))

generate_half_zeroes_half_ones("../meshes/half3d.mesh", (32, 32, 32), (10.0, 10.0, 10.0))
generate_half_zeroes_half_ones("../meshes/half2d.mesh", (64, 64), (10.0, 10.0))

generate_box_2d("../meshes/box2d_512.mesh", (512, 512), (10.0, 10.0), 4.0)
generate_box_3d("../meshes/box3d.mesh", (32, 32, 32), (10.0, 10.0, 10.0), 4.0)
