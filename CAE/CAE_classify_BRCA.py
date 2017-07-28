from keras.models import Sequential
from keras.layers.core import Dense, Dropout, Activation, Flatten, Reshape
from keras.layers.convolutional import Convolution2D, MaxPooling2D, UpSampling2D
from keras.layers.normalization import BatchNormalization
from keras.layers import Dropout
from keras.optimizers import SGD, Adam
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
import matplotlib
import pickle
from matplotlib import pyplot as plt
from matplotlib import offsetbox
import numpy as np
import os.path
import image
import copy
import PIL
import sys
import PIL.ImageOps
import ipdb

import warnings
warnings.filterwarnings("ignore")

weightsFile = 'CAE_BRCA/pretrained_7999.h5'
classFile = 'BRCA_classes.txt'
imDir = 'BRCA_images/'
newWeightsFile = 'CAE_classify_BRCA/discrim_weights'
lossFile = 'CAE_classify_BRCA_loss.txt'

batchSize = 1000 # number of samples to train per epoch
numEpochs = 3000 # number of network training iterations
newSize   = 512 # resize images before sampling
sampSize  = 128 # size of sampling window after resizing
filterSize= 5 # conv kernel size
poolSize  = 2 # pooling size

trnFrac = float(sys.argv[1]) # fraction for training
valFrac = float(sys.argv[2]) # fraction for validation
trnBalance = '1' == sys.argv[3] # balance training set classes?

print('-' * 50)

print('Training Fraction:   ' + str(trnFrac))
print('Validation Fraction: ' + str(valFrac))
print('Testing Fraction:    ' + str(1 - trnFrac - valFrac))
print('-' * 50)
if trnBalance:
	print('Balancing Training Set')
else:
	print('Not Balancing Training Set')
print('-' * 50)

allPaths = np.asarray(image.list_pictures(imDir))
np.random.shuffle(allPaths)
numIms = len(allPaths)
trainIms = allPaths[:numIms * trnFrac]
testIms = allPaths[numIms * (trnFrac + valFrac):]
valIms = allPaths[numIms * trnFrac : numIms * (trnFrac + valFrac)]

# pull sample classes
labels = {}
classNames = {}
classIms = {}
classCount = 0
fi = open(classFile)
for lin in fi:
	im, cl = lin.strip('\n').split('\t')
	if cl in classNames:
		classNum = classNames[cl]
	else:
		classNames[cl] = classCount
		classNum = classCount
		classIms[classCount] = []
		classCount += 1
	classIms[classNum].append(im)
	labels[im] =  classNum

trnClassIms = {}
for c in xrange(classCount):
	trnClassIms[c] = []

for im in labels:
	if (imDir + im) in trainIms:
		trnClassIms[labels[im]].append(im)

nClasses = len(classNames)
for c in classIms:
	classIms[c] = np.asarray(classIms[c])
	trnClassIms[c] = np.asarray(trnClassIms[c])

sums = np.zeros((nClasses))
for c in xrange(nClasses):
	sums[c] = len(trnClassIms[c])

if sum(np.asarray(sums) >= 0) != nClasses:
	print('not all classes in training set. rerun.')
	sys.exit()

indsLow = np.argsort(sums)
toAdd = 0
for i in xrange(nClasses - 1):
	toAdd += max(sums) - sums[indsLow[i]]

imsAdd = np.asarray([0] * toAdd).astype('S64')  # changed 'string' to 'S64' -DM
soFar = 0
for i in xrange(nClasses - 1):
	amount = max(sums) - sums[indsLow[i]]
	inds = np.floor(np.random.rand((amount)) * sums[indsLow[i]]).astype('int')
	imsAdd[soFar:soFar+amount] = trnClassIms[indsLow[i]][inds]
	soFar += amount

co = len(trainIms)
tmp = np.zeros((co + toAdd)).astype('S64')  # changed 'string' to 'S64' -DM
tmp[:co] = trainIms
for im in imsAdd:
	tmp[co] = imDir + im
	co += 1

if trnBalance:
	trainIms = tmp

# helper function to get a random sampling from an image directory
def getBatch(batchSize, allPaths, numClass):
	allIms = np.zeros((batchSize, 3, sampSize, sampSize))
	Y = np.zeros((batchSize, 1)).astype('int')
	rIms = np.zeros((batchSize, 1)).astype('str')
	for i in xrange(batchSize):
		f1 = int(np.sign(np.random.rand() - .5))
		f2 = int(np.sign(np.random.rand() - .5))
		rIm = allPaths[int(np.random.rand() * len(allPaths))]
		rIms[i] = rIm
		im = image.load_img(rIm).resize((newSize, newSize))
		Y[i] = labels[rIm.split('/')[-1]]
		r  = int(np.random.rand() * (newSize - sampSize))
		c  = int(np.random.rand() * (newSize - sampSize))
		im = im.crop((r, c, r + sampSize, c + sampSize))
		im = image.random_rotation(im, 5)
		allIms[i, :, :, :] = image.img_to_array(im)
		allIms[i, :, :, :] = allIms[i, :, ::-f1, ::-f2]
	if numClass > 2:
		Y2 = np.zeros((batchSize, numClass)).astype('int')
		for i in xrange(batchSize):
			Y2[i, Y[i][0]] = 1
		Y = Y2
	return allIms, Y, rIms

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
cae.add(Dense(128 * 4 * 4))
cae.add(Activation('relu'))
cae.add(Reshape((128, 4, 4)))
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
cae.load_weights(weightsFile)

# initialize classifier with pretrained weights
cl = Sequential()
cl.add(Convolution2D(8, filterSize, filterSize,
	input_shape=(3, sampSize, sampSize), border_mode='same',  weights=cae.layers[0].get_weights()))
cl.add(MaxPooling2D(pool_size=(poolSize, poolSize)))
cl.add(Activation('relu'))
cl.add(Convolution2D(16, filterSize, filterSize, border_mode='same', weights=cae.layers[3].get_weights()))
cl.add(MaxPooling2D(pool_size=(poolSize, poolSize)))
cl.add(Activation('relu'))
cl.add(Convolution2D(32, filterSize, filterSize, border_mode='same', weights=cae.layers[6].get_weights()))
cl.add(MaxPooling2D(pool_size=(poolSize, poolSize)))
cl.add(Activation('relu'))
cl.add(Convolution2D(64, filterSize, filterSize, border_mode='same', weights=cae.layers[9].get_weights()))
cl.add(MaxPooling2D(pool_size=(poolSize, poolSize)))
cl.add(Activation('relu'))
cl.add(Convolution2D(128, filterSize, filterSize, border_mode='same', weights=cae.layers[12].get_weights()))
cl.add(MaxPooling2D(pool_size=(poolSize, poolSize)))
cl.add(Activation('relu'))
cl.add(Flatten())
cl.add(Dense(1024, weights=cae.layers[16].get_weights()))
cl.add(Activation('relu'))
sgd = SGD(lr=0.1, decay=1e-6, momentum=0.9, nesterov=True)


cl.add(Dense(100))
cl.add(Dropout(0.5))
cl.add(Activation('relu'))

if nClasses == 2:
	cl.add(Dense(1))
	cl.add(Dropout(0.5))
	cl.add(Activation('sigmoid'))
	cl.compile(loss='binary_crossentropy', optimizer=Adam(lr=5e-7, decay=1e-5))
else:
	cl.add(Dense(nClasses))
	cl.add(Dropout(0.5))
	cl.add(Activation('softmax'))
	cl.compile(loss='categorical_crossentropy', optimizer=Adam(lr=5e-7, decay=1e-5))


for i in xrange(numEpochs):
	X, Y, tmp = getBatch(batchSize, trainIms, nClasses)
	# Xval, Yval, tmp = getBatch(batchSize, valIms, nClasses)
	# cl.fit(X, Y, validation_data=(Xval, Yval), nb_epoch=1, verbose=0)
	cl.fit(X, Y, nb_epoch=1, verbose=0)
	loss= cl.evaluate(X, Y, verbose=0)
	open(lossFile, 'a').write(str(loss) + '\n')
	if (i + 1) % 50 == 0:
                cl.save_weights(newWeightsFile + '_' + str(i) + '.h5')
