import datetime

class sequence:
	def __init__(self):
		self.l = [1,1,1]
		self.fitness = 0
	def dump(self,path):
		path = path + '/{:%Y%m%d%H%M%S}.tsv'.format(datetime.datetime.now())
		with open(path,'w+') as f:
			f.write('hello')
