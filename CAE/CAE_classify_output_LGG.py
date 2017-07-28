from keras.models import Sequential
from keras.layers.core import Dense, Dropout, Activation, Flatten, Reshape
from keras.layers.convolutional import Convolution2D, MaxPooling2D, UpSampling2D
from keras.layers.normalization import BatchNormalization
from sklearn.decomposition import PCA
import numpy as np
import os.path
import image
import copy
import PIL
import PIL.ImageOps

weightsFile = 'CAE_classify_LGG/discrim_weights_2999.h5' # trained network weights
outputFile  = 'CAE_classify_output_LGG.csv' # name of file to save
imPath      = 'LGG_images/' # path to images
nClasses = 2

batchSize   = 1000 # number of samples to train per epoch
numEpochs   = 3000 # number of network training iterations
newSize     = 512 # resize images before sampling
sampSize    = 128 # size of sampling window after resizing

filterSize  = 5 # conv kernel size
poolSize    = 2 # pooling size

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
		allIms[i, :, :, :] = allIms[i, :, ::-f1, ::-f2] / 255.
	return allIms

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
		allIms[i, :, :, :] = allIms[i, :, ::-f1, ::-f2] / 255.
	return allIms

def PIL2array(img):
	return np.array(img.getdata()).reshape(img.size[1], img.size[0], 3)

# initialize classifier with pretrained weights
cl = Sequential()
cl.add(Convolution2D(8, filterSize, filterSize,
	input_shape=(3, sampSize, sampSize), border_mode='same'))
cl.add(MaxPooling2D(pool_size=(poolSize, poolSize)))
cl.add(Activation('relu'))
cl.add(Convolution2D(16, filterSize, filterSize, border_mode='same'))
cl.add(MaxPooling2D(pool_size=(poolSize, poolSize)))
cl.add(Activation('relu'))
cl.add(Convolution2D(32, filterSize, filterSize, border_mode='same'))
cl.add(MaxPooling2D(pool_size=(poolSize, poolSize)))
cl.add(Activation('relu'))
cl.add(Convolution2D(64, filterSize, filterSize, border_mode='same'))
cl.add(MaxPooling2D(pool_size=(poolSize, poolSize)))
cl.add(Activation('relu'))
cl.add(Convolution2D(128, filterSize, filterSize, border_mode='same'))
cl.add(MaxPooling2D(pool_size=(poolSize, poolSize)))
cl.add(Activation('relu'))
cl.add(Flatten())
cl.add(Dense(1024))
cl.add(Activation('relu'))
cl.add(Dense(100))
cl.add(Dropout(0.5))
cl.add(Activation('relu'))

if nClasses == 2:
	cl.add(Dense(1))
	cl.add(Dropout(0.5))
	cl.add(Activation('sigmoid'))
	cl.compile(loss='binary_crossentropy', optimizer='adam')
else:
	cl.add(Dense(nClasses))
	cl.add(Dropout(0.5))
	cl.add(Activation('softmax'))
	cl.compile(loss='categorical_crossentropy', optimizer='adam')

cl.load_weights(weightsFile)

encode = Sequential()
encode.add(Convolution2D(8, filterSize, filterSize,
	input_shape=(3, sampSize, sampSize), border_mode='same',  weights=cl.layers[0].get_weights()))
encode.add(MaxPooling2D(pool_size=(poolSize, poolSize)))
encode.add(Activation('relu'))
encode.add(Convolution2D(16, filterSize, filterSize, border_mode='same', weights=cl.layers[3].get_weights()))
encode.add(MaxPooling2D(pool_size=(poolSize, poolSize)))
encode.add(Activation('relu'))
encode.add(Convolution2D(32, filterSize, filterSize, border_mode='same', weights=cl.layers[6].get_weights()))
encode.add(MaxPooling2D(pool_size=(poolSize, poolSize)))
encode.add(Activation('relu'))
encode.add(Convolution2D(64, filterSize, filterSize, border_mode='same', weights=cl.layers[9].get_weights()))
encode.add(MaxPooling2D(pool_size=(poolSize, poolSize)))
encode.add(Activation('relu'))
encode.add(Convolution2D(128, filterSize, filterSize, border_mode='same', weights=cl.layers[12].get_weights()))
encode.add(MaxPooling2D(pool_size=(poolSize, poolSize)))
encode.add(Activation('relu'))
encode.add(Flatten())
encode.add(Dense(1024, weights=cl.layers[16].get_weights()))
encode.add(Activation('relu'))
encode.add(Dense(100, weights=cl.layers[18].get_weights()))
encode.add(Dropout(0.5))
encode.add(Activation('relu'))
encode.compile(loss='mse', optimizer='adam')

###############################

allPaths = image.list_pictures(imPath)
repAll = np.zeros((len(allPaths), 100))

for i in xrange(len(repAll)):
	for co in xrange(10):
		repAll[i, :] += sum(encode.predict(getOne(100, allPaths[i])))
	print(i)
repAll /= 1000.

K = 100
pca = PCA(n_components=K, whiten=True)
reduced = pca.fit_transform(repAll)
vari = pca.explained_variance_ratio_
nDims = np.where(np.cumsum(vari) > .98)[0][0]
# repAll = reduced[:,:nDims]
nDims = reduced.shape[1]

fi = open(outputFile,'w')
for i in xrange(len(repAll)):
        fi.write(os.path.basename(allPaths[i]) + ',')
        for dim in xrange(nDims):
                fi.write(str(repAll[i][dim]))
                if dim < nDims - 1:
                        fi.write(',')
        fi.write('\n')
fi.close()
