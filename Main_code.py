# -*- coding: utf-8 -*-
"""Untitled23.ipynb

Automatically generated by Colab.

Original file is located at
    https://colab.research.google.com/drive/1HSfljRHzpyLzWsKqTvpsf2cchy4UqAYC
"""

################# installing of libraries ############################

!pip install tensorflow
!pip install torch torchvision scikit-learn
!pip install imbalanced-learn


############# load packages ##########################################

import numpy as np
import pandas as pd
import statistics
import matplotlib.pyplot as plt
from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import StandardScaler
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from sklearn.model_selection import KFold
from sklearn import metrics
from sklearn.metrics import precision_score,f1_score,average_precision_score,recall_score
from sklearn.metrics import precision_recall_fscore_support
from imblearn.metrics import specificity_score
from keras.models import Sequential
from keras.layers import LSTM
from keras.layers import Dense
from keras.layers import RepeatVector
from keras.layers import TimeDistributed



############## encoding of features and class attributes ##############

# Loading Data
sample_df=pd.read_csv('File_name.csv')
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


############## Stacked LSTM Autoencoder Model #############

inputs = keras.layers.Input(shape=(X.shape[1],1))
input_dim=X.shape[1]
losses_autoencoder=[]

# Define encoder

encoder = LSTM(150, activation='relu',return_sequences=True)(inputs)
encoder=LSTM(50, activation='relu', return_sequences=False)(encoder)
encoder1 = RepeatVector(X.shape[1])(encoder)

# Define decoder

decoder1 = LSTM(50, activation='relu', return_sequences=True)(encoder1)
decoder1 = LSTM(150, activation='relu', return_sequences=True)(decoder1)
output = TimeDistributed(Dense(1))(decoder1)

# Tie it together
autoencoder = keras.Model(inputs, output)

# Compilation of model
autoencoder.compile(optimizer=tf.keras.optimizers.Adam(learning_rate=0.001), loss='mse')

# Fit of model
history_autoencoder = autoencoder.fit(X, X,batch_size=200,epochs=40,validation_split=0.2, verbose=1)


# Measuring loss
losses_autoencoder.append(history_autoencoder.history['loss'])
print(losses_autoencoder)

# Execution of encoder to fetch important features
encoder2 = keras.Model(inputs,encoder)
yhat = encoder2.predict(X_feat_scaled)


########## defining classifier model ##################################

num_classes = len(np.unique(y))

# Define the number of folds for cross validation
num_folds = 10

# Initialize variables to store overall performance metrics
accuracies = []
losses_autoencoder = []
losses_classification = []
precision_score=[]
recall_score=[]
f1score=[]
sp_score=[]


# Initialize the KFold object
kf = KFold(n_splits=num_folds, shuffle=True, random_state=42)

for train_index, test_index in kf.split(yhat):
    X_train, X_test = yhat[train_index], yhat[test_index]
    y_train, y_test = y[train_index], y[test_index]

    classification_model = keras.Sequential([
        layers.Dense(35, activation="relu", input_shape=(yhat.shape[1],)),
        layers.Dense(25, activation="relu", input_shape=(yhat.shape[1],)),
        layers.Dense(num_classes, activation="softmax")
    ])

    classification_model.compile(optimizer="adam", loss="sparse_categorical_crossentropy", metrics=["accuracy"])


    # Train the classification model and get the training history
    history_classification = classification_model.fit(X_train, y_train,batch_size=64,epochs=200, validation_split=0.2, verbose=1)
    losses_classification.append(history_classification.history['loss'])

    # Evaluation of the classifier model
    accuracy = classification_model.evaluate(X_test, y_test, verbose=0)[1]

    # Evaluation of Accuracy
    accuracies.append(accuracy)

    # Evaluation of Precision,Recall,F1_score
    y_pred = tf.argmax(classification_model.predict(X_test), axis=1)

    precision_metric = metrics.precision_score(y_test,y_pred,labels=None, pos_label=1, average ='weighted',sample_weight=None, zero_division='warn')
    recall_metric=metrics.recall_score(y_test, y_pred,labels=None, pos_label=1, average='weighted', sample_weight=None, zero_division='warn')
    f1_metric=metrics.f1_score(y_test, y_pred,labels=None, pos_label=1, average='weighted', sample_weight=None, zero_division='warn')
    sp_metric=specificity_score(y_test, y_pred, average='weighted')

    precision_score.append(precision_metric)
    recall_score.append(recall_metric)
    f1score.append(f1_metric)
    sp_score.append(sp_metric)

# Printing average of overall performance matrics
avg_p=np.mean(precision_score)
avg_r=np.mean(recall_score)
avg_f=np.mean(f1score)
avg_sp=np.mean(sp_score)
average_accuracy = np.mean(accuracies)

print("Average Accuracy: {:.5f}%".format(average_accuracy * 100))
print("Average Precision:{:.5f}%".format(avg_p))
print("Average Recall:{:.5f}%".format(avg_r))
print("Average F1_score:{:.5f}%".format(avg_f))
print("Average Specificity:{:.5f}%".format(avg_sp))
print("Deviation of Precision Value",statistics.stdev(precision_score))
print("Deviation of Recall Value",statistics.stdev(recall_score))
print("Deviation of F1Score Value",statistics.stdev(f1score))
print("Deviation of Specificity Value",statistics.stdev(sp_score))



############ validation on new set of data ################################


# Reading of data
sample_v = pd.read_csv('validation_data.csv')

X_originalv = sample_v.drop('Label', axis=1)
y_originalv = sample_v['Label']

features_v = pd.DataFrame(X_originalv)
labels_v = pd.DataFrame(y_originalv)
X_featv = features_v.values


# Apply standard scaler to all features
scaler = StandardScaler()
X_feat_scaledv = scaler.fit_transform(X_featv)

# Apply label encoder to all labels
y_labelsv=labels_v.values


# Final X (feature) and y (class-label) matrices for validation data

X_v=X_feat_scaledv
yLv=y_labelsv.ravel()
y_v=yLv

# Fetching of important features
yhatv = encoder2.predict(X_v)

# Evaluation of Accuracy
y_predv = tf.argmax(classification_model.predict(yhatv), axis=1)
val_accuracy = classification_model.evaluate(yhatv, y_v, verbose=0)[1]

print("Accuracy",val_accuracy)

# Evaluation of performance matrix
precision_metricv = metrics.precision_score(y_v,y_predv,labels=None, pos_label=1, average ='weighted',sample_weight=None, zero_division='warn')
recall_metricv=metrics.recall_score(y_v, y_predv,labels=None, pos_label=1, average='weighted', sample_weight=None, zero_division='warn')
f1_sv=metrics.f1_score(y_v, y_predv,labels=None, pos_label=1, average='weighted')
sp_v=specificity_score(y_v, y_predv, average='weighted')

print("Precision of validated Data",precision_metricv)
print("Recall of validated Data",recall_metricv)
print("F1_Score of validated Data",f1_sv)
print("Specificity of validated Data",sp_v)