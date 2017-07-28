from keras.models import Sequential
from keras.layers.core import Dense, Dropout, Activation, Flatten, Reshape
from keras.layers.convolutional import Convolution2D, MaxPooling2D, UpSampling2D
from keras.layers.normalization import BatchNormalization
from keras.optimizers import Adam
import numpy as np
import os.path
import image

batchSize = 1000 # number of samples to train per epoch
startEpoch = 0
numEpochs = 8000 # number of network training iterations
newSize  = 512 # resize images before sampling
sampSize = 128 # size of sampling window after resizing

filterSize = 5 # conv kernel size
poolSize   = 2 # pooling size

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
	return allIms

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
cae.add(Activation('sigmoid'))  # ADDITION -DM

# compile and load pretrained weights
cae.compile(loss='mse', optimizer=Adam(lr=0.0005, decay=1e-5))
if os.path.isfile('CAE_BRCA/pretrained.h5'):
	cae.load_weights('CAE_BRCA/pretrained.h5')

# train the model
for i in xrange(startEpoch, startEpoch + numEpochs):
	imBatch = getBatch(batchSize, 'BRCA_images/') / 255.
	cae.fit(imBatch, imBatch, nb_epoch=1, verbose=True, batch_size=5)
#	cae.save_weights('pretrained.h5',overwrite=True)
	# also save the weights every 20 iterations just in case training gets cut off
	if (i + 1) % 50 == 0:
		cae.save_weights('CAE_BRCA/pretrained_' + str(i) + '.h5')
