from keras.models import Sequential
from keras.layers.core import Dense, Dropout, Activation, Flatten, Reshape
from keras.layers.convolutional import Convolution2D, MaxPooling2D, UpSampling2D
from keras.layers.normalization import BatchNormalization
from sklearn.decomposition import PCA
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import offsetbox
import numpy as np
import os.path
import image
import copy
import PIL
import PIL.ImageOps

imagePath = 'GTEx_images/'
pretrained = 'CAE_GTEx/pretrained_11999.h5'
outputName = 'CAE_output_GTEx' # create this folder before running

batchSize = 1000 # number of samples to train per epoch
numEpochs = 3000 # number of network training iterations
newSize  = 512 # resize images before sampling
sampSize = 128 # size of sampling window after resizing

filterSize = 5 # conv kernel size
poolSize = 2 # pooling size

# helper function to get a random sampling from an image directory
def getBatch(batchSize, path):
	allPaths = image.list_pictures(path)
	allIms = np.zeros((batchSize, 3, sampSize, sampSize))
	for i in xrange(batchSize):
		f1 = int(np.sign(np.random.rand() - .5))
		f2 = int(np.sign(np.random.rand() - .5))
		im = image.load_img(allPaths[int(np.random.rand() * len(allPaths))]).resize((newSize, newSize))
		r  = int(np.random.rand() * (newSize - sampSize))
		c  = int(np.random.rand() * (newSize - sampSize))
		im = im.crop((r, c, r + sampSize, c + sampSize))
		im = image.random_rotation(im, 5)
		allIms[i, :, :, :] = image.img_to_array(im)
		allIms[i, :, :, :] = allIms[i, :, ::-f1, ::-f2]
	return allIms/255.

def getOne(batchSize, imPath):
	allIms = np.zeros((batchSize, 3, sampSize, sampSize))
	for i in xrange(batchSize):
		f1 = int(np.sign(np.random.rand() - .5))
		f2 = int(np.sign(np.random.rand() - .5))
		im = image.load_img(imPath).resize((newSize, newSize))
		r  = int(np.random.rand() * (newSize - sampSize))
		c  = int(np.random.rand() * (newSize - sampSize))
		im = im.crop((r, c, r + sampSize, c + sampSize))
		im = image.random_rotation(im, 5)
		allIms[i, :, :, :] = image.img_to_array(im)
		allIms[i, :, :, :] = allIms[i, :, ::-f1, ::-f2]
	return allIms/255.

def plot_embedding(X, ims, title=None):
	x_min, x_max =  np.min(X, 0), np.max(X, 0)
	X = (X - x_min) / (x_max - x_min)
	plt.figure()
	ax = plt.subplot(111)
	for i in range(X.shape[0]):
		plt.plot(X[i, 0], X[i, 1], 'b.')
	if hasattr(offsetbox, 'AnnotationBbox'):
        # only print thumbnails with matplotlib > 1.0
		shown_images = np.array([[1., 1.]])  # just something big
		for i in range(ims.shape[0]):
			dist = np.sum((X[i] - shown_images) ** 2, 1)
			if np.min(dist) < 4e-3:
		                # don't show points that are too close
				continue
			shown_images = np.r_[shown_images, [X[i]]]
			imagebox = offsetbox.AnnotationBbox(
				offsetbox.OffsetImage(ims[i]),
				X[i],pad=0)
			ax.add_artist(imagebox)
	plt.xticks([]), plt.yticks([])
	if title is not None:
		plt.title(title)

def PIL2array(img):
	return np.array(img.getdata()).reshape(img.size[1], img.size[0], 3)

# initialize cae
cae = Sequential()

# convolution + pooling 1
cae.add(Convolution2D(8, filterSize, filterSize, input_shape=(3, sampSize, sampSize), border_mode='same'))
cae.add(MaxPooling2D(pool_size=(poolSize, poolSize)))
cae.add(Activation('relu'))

# convolution + pooling 2
cae.add(Convolution2D(16, filterSize, filterSize, border_mode='same'))
cae.add(MaxPooling2D(pool_size=(poolSize, poolSize)))
cae.add(Activation('relu'))

# convolution + pooling 3
cae.add(Convolution2D(32, filterSize, filterSize, border_mode='same'))
cae.add(MaxPooling2D(pool_size=(poolSize, poolSize)))
cae.add(Activation('relu'))

# convolution + pooling 4
cae.add(Convolution2D(64, filterSize, filterSize, border_mode='same'))
cae.add(MaxPooling2D(pool_size=(poolSize, poolSize)))
cae.add(Activation('relu'))

# convolution + pooling 5
cae.add(Convolution2D(128, filterSize, filterSize, border_mode='same'))
cae.add(MaxPooling2D(pool_size=(poolSize, poolSize)))
cae.add(Activation('relu'))

# dense network
cae.add(Flatten())
cae.add(Dense(1024))
cae.add(Activation('relu'))
cae.add(Dense(128*4*4))
cae.add(Activation('relu'))
cae.add(Reshape(dims=(128, 4, 4)))
cae.add(Activation('relu'))

# unpooling + deconvolution 1
cae.add(UpSampling2D(size=(poolSize, poolSize)))
cae.add(Convolution2D(64, filterSize, filterSize, border_mode='same'))
cae.add(Activation('relu'))

# unpooling + deconvolution 2
cae.add(UpSampling2D(size=(poolSize, poolSize)))
cae.add(Convolution2D(32, filterSize, filterSize, border_mode='same'))
cae.add(Activation('relu'))

# unpooling + deconvolution 3
cae.add(UpSampling2D(size=(poolSize, poolSize)))
cae.add(Convolution2D(16, filterSize, filterSize, border_mode='same'))
cae.add(Activation('relu'))

# unpooling + deconvolution 4
cae.add(UpSampling2D(size=(poolSize, poolSize)))
cae.add(Convolution2D(8, filterSize, filterSize, border_mode='same'))
cae.add(Activation('relu'))

# final unpooling + deconvolution
cae.add(UpSampling2D(size=(poolSize, poolSize)))
cae.add(Convolution2D(3, filterSize, filterSize,  border_mode='same'))
cae.add(Activation('sigmoid'))

# compile and load pretrained weights
cae.compile(loss='mse', optimizer='adam')
cae.load_weights(pretrained)

# initialize encoder
encode = Sequential()
encode.add(Convolution2D(8, filterSize, filterSize,
	input_shape=(3, sampSize, sampSize), border_mode='same',  weights=cae.layers[0].get_weights()))
encode.add(MaxPooling2D(pool_size=(poolSize, poolSize)))
encode.add(Activation('relu'))
encode.add(Convolution2D(16, filterSize, filterSize, border_mode='same', weights=cae.layers[3].get_weights()))
encode.add(MaxPooling2D(pool_size=(poolSize, poolSize)))
encode.add(Activation('relu'))
encode.add(Convolution2D(32, filterSize, filterSize, border_mode='same', weights=cae.layers[6].get_weights()))
encode.add(MaxPooling2D(pool_size=(poolSize, poolSize)))
encode.add(Activation('relu'))
encode.add(Convolution2D(64, filterSize, filterSize, border_mode='same', weights=cae.layers[9].get_weights()))
encode.add(MaxPooling2D(pool_size=(poolSize, poolSize)))
encode.add(Activation('relu'))
encode.add(Convolution2D(128, filterSize, filterSize, border_mode='same', weights=cae.layers[12].get_weights()))
encode.add(MaxPooling2D(pool_size=(poolSize, poolSize)))
encode.add(Activation('relu'))
encode.add(Flatten())
encode.add(Dense(1024, weights=cae.layers[16].get_weights()))
encode.add(Activation('relu'))
encode.compile(loss='mse', optimizer='adam')

# initialize decoder
decode = Sequential()
decode.add(Dense(128*4*4, input_dim=(1024), weights=cae.layers[18].get_weights()))
decode.add(Activation('relu'))
decode.add(Reshape(dims=(128, 4, 4)))
decode.add(Activation('relu'))
decode.add(UpSampling2D(size=(poolSize, poolSize)))
decode.add(Convolution2D(64, filterSize, filterSize, border_mode='same', weights=cae.layers[23].get_weights()))
decode.add(Activation('relu'))
decode.add(UpSampling2D(size=(poolSize, poolSize)))
decode.add(Convolution2D(32, filterSize, filterSize, border_mode='same', weights=cae.layers[26].get_weights()))
decode.add(Activation('relu'))
decode.add(UpSampling2D(size=(poolSize, poolSize)))
decode.add(Convolution2D(16, filterSize, filterSize, border_mode='same', weights=cae.layers[29].get_weights()))
decode.add(Activation('relu'))
decode.add(UpSampling2D(size=(poolSize, poolSize)))
decode.add(Convolution2D(8, filterSize, filterSize, border_mode='same', weights=cae.layers[32].get_weights()))
decode.add(Activation('relu'))
decode.add(UpSampling2D(size=(poolSize, poolSize)))
decode.add(Convolution2D(3, filterSize, filterSize, border_mode='same', weights=cae.layers[35].get_weights()))
decode.add(Activation('sigmoid'))
decode.compile(loss='mse', optimizer='adam')



# #############################


allPaths = image.list_pictures(imagePath)

repAll = np.zeros((len(allPaths), 1024))
for i in xrange(len(repAll)):
	print(i)
	for co in xrange(10):
		repAll[i, :] += sum(encode.predict(getOne(100, allPaths[i])))
	repAll[i, :] /= 1000.


# this shows the encodeing and decoding for some random images
ims = getBatch(101, imagePath)
for i in xrange(100):
	image.array_to_img(ims[i] * 255.).save(outputName + '/' + str(i) + 'A.png')
	image.array_to_img(cae.predict(ims[i:i+1])[0] * 255.).save(outputName + '/' + str(i) + 'B.png')

K = 1024
pca = PCA(n_components=K, whiten=True)
reduced = pca.fit_transform(repAll)
vari = pca.explained_variance_ratio_
# nDims = np.where(np.cumsum(vari) > .98)[0][0]
# reduced = reduced[:,:nDims]

nDims = reduced.shape[1]

with open(outputName + '_vari.txt', 'w') as fi:
        for x in np.cumsum(vari):
                fi.write(str(x) + '\n')

fi = open(outputName + '_rep.csv','w')
for i in xrange(len(repAll)):
        fi.write(os.path.basename(allPaths[i]) + ',')
        for dim in xrange(nDims):
                fi.write(str(repAll[i][dim]))
                if dim < nDims - 1:
                        fi.write(',')
        fi.write('\n')
fi.close()

fi = open(outputName + '_repPCA.csv','w')
for i in xrange(len(repAll)):
        fi.write(os.path.basename(allPaths[i]) + ',')
	for dim in xrange(nDims):
		fi.write(str(reduced[i][dim]))
		if dim < nDims - 1:
			fi.write(',')
	fi.write('\n')
fi.close()
