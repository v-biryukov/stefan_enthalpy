import struct


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




generate_zeroes("../meshes/zeroes3d.mesh", (32, 32, 32), (10.0, 10.0, 10.0))
generate_zeroes("../meshes/zeroes2d.mesh", (64, 64), (10.0, 10.0))

generate_half_zeroes_half_ones("../meshes/half3d.mesh", (32, 32, 32), (10.0, 10.0, 10.0))
generate_half_zeroes_half_ones("../meshes/half2d.mesh", (64, 64), (10.0, 10.0))