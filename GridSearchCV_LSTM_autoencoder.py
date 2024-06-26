# -*- coding: utf-8 -*-
"""Untitled24.ipynb

Automatically generated by Colab.

Original file is located at
    https://colab.research.google.com/drive/1KzvblVZNhU59FrE4zICGjJyLDmVxasx-
"""

################# installing of libraries ######################################

!pip install tensorflow
!pip install torch torchvision scikit-learn
!pip install scikeras


################ load packages #################################################

import numpy as np
import pandas as pd
import tensorflow as tf
import matplotlib.pyplot as plt
from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import StandardScaler
from tensorflow import keras
from tensorflow.keras import layers
from keras.models import Sequential
from keras.layers import LSTM
from keras.layers import Dense
from keras.layers import RepeatVector
from keras.layers import TimeDistributed
from sklearn import preprocessing
from sklearn.model_selection import GridSearchCV
from scikeras.wrappers import KerasClassifier


######################### input data and pre-processing  #########################

# Loading Data
sample_df=pd.read_csv('Acinetobacter_4mer.csv')
X_original = sample_df.drop('Label', axis=1)
y_original = sample_df['Label']

# Encoding of features and class Labels
features = pd.DataFrame(X_original)
labels = pd.DataFrame(y_original)
X_feat = features.values

# Apply standard scaler to all features
scaler = StandardScaler()
X_feat_scaled = scaler.fit_transform(X_feat)

# Apply label encoder to all labels
y_labels=labels.values
label_encoder = preprocessing.LabelEncoder()
yL= label_encoder.fit_transform(y_labels)

# Final X (feature) and y (class-label) matrices
y = yL
X = X_feat_scaled


#######################  defining Stacked LSTM Autoencoder model and applying GridSearchCV ###############

inputs = keras.layers.Input(shape=(X.shape[1],1))
input_dim=X.shape[1]
losses_autoencoder=[]

def create_model(neurons1,neurons2,lr):
# Defining encoder and decoder
    model = Sequential()
    model.add(LSTM(neurons1, activation='relu', input_shape=(X.shape[1],1),return_sequences=True))
    model.add(LSTM(neurons2, activation='relu',return_sequences=False))
    model.add(RepeatVector(X.shape[1]))
    model.add(LSTM(neurons2, activation='relu', return_sequences=True))
    model.add(LSTM(neurons1, activation='relu',return_sequences=True))
    model.add(TimeDistributed(Dense(1)))
    model.compile(optimizer=tf.keras.optimizers.Adam(learning_rate=lr), loss='mse')
    return model


# Fitting the model and getting the best estimator
model=KerasClassifier(build_fn=create_model,neurons1=[50,100,150,200],neurons2=[20,30,40,50,60],lr=[0.1,0.01,0.001])

params={ 'lr':[0.1,0.01,0.001],
         'neurons1':[50,100,150,200],
         'neurons2':[20,30,40,50,60]
        }

gs=GridSearchCV(estimator=model, param_grid=params, cv=5,scoring='accuracy')

# now fit the dataset to the GridSearchCV object.
gs_result = gs.fit(X, y)
print("Best: %f using %s" % (gs_result.best_score_, gs_result.best_params_))