#!/usr/bin/Rscript
#Cours: Deep Learning
#Prof: Frederic Guyon
#Octobre 2018


#Université Paris Diderot ~ Master Biologie Informatique M2


#import libraries
library(keras)

# 
# # génération des données à partir modèle y=cos(x)+0.1*x^2
# x <- seq(-5,5,len=50)
# ytrue <- cos(x)+0.1*x*x
# yobs <- ytrue+rnorm(length(x),sd = 0.2)
# #visualisation
# plot(x,ytrue,ty='l',col="red") #real values
# points(x,yobs, col="blue", pch="+") #variations
# 
# #Neural Network
# model = keras_model_sequential()
# layer_dense(model, units=20, activation="relu", input_shape=20)
# layer_dense(model, units=1, activation="linear", input_shape=20)
# compile(model,loss = "mean_squared_error", optimizer = optimizer_rmsprop())
# fit(model, x, yobs, epochs = 1500, batch_size = 10)
# evaluate(model, x, yobs)
# ypred=predict(model,x)
# points(x,ypred, col="red")
# 
# 
# #*********************************************************
# #Exerice_1 ~ Classification de données chiffres manuscrit
# #********************************************************
# #1/
# mnist <- dataset_mnist()
# 
# 
# #2/visualise some minst data
# par(mfrow = c(2,2))
# for (i in 1:4){
#   img <- t(apply(mnist$train$x[i,,], 2, rev))
#   image(img,col=gray(0:255/255))
#   print(mnist$train$y[i])
# }
# dev.off()
# 
# #3/Apprentissage et Validation des données
# #Pretreament data
# x_train <- mnist$train$x
# y_train <- mnist$train$y
# x_test <- mnist$test$x
# y_test <- mnist$test$y
# # reshape
# x_train <- array_reshape(x_train, c(nrow(x_train), 784))
# x_test <- array_reshape(x_test, c(nrow(x_test), 784))
# # rescale
# x_train <- x_train / 255
# x_test <- x_test / 255
# # output generation
# y_train <- to_categorical(y_train, 10)
# y_test <- to_categorical(y_test, 10)
# 
# #Processing
# model = keras_model_sequential()
# 
# #Training
# #input_shape obligatoire comme premiere couche, dim. output=256
# model = layer_dense(model,units = 256, activation = "relu", input_shape = 784) #units ~ neural layer
# #input_shape 
# # model = layer_dropout(model,rate = 0.4) 
# # model = layer_dense(model,units = 128, activation = "relu")
# model = layer_dropout(model,rate = 0.3)
# model = layer_dense(model,units = 10, activation = "softmax")
# model=compile(model, loss = "categorical_crossentropy", optimizer = optimizer_rmsprop(), metrics = c("accuracy"))
# history=fit(model, x_train, y_train, epochs = 30, batch_size = 128, validation_split = 0.2)
# #visualisation
# plot(history$metrics$acc)
# points(history$metrics$val_acc,col="blue")
# 
# #Validation
# y_pred=predict_classes(model, x_test)
# 
# #Evaluation
# evaluate(model, x_test, y_test)
# Table=table(y_pred, mnist$test$y)
# sum(diag(Table))/10000


#*************************************************
#Exercice 2 : Classification des données Fashions
#*************************************************
#loading dataset
fashion_mnist = dataset_fashion_mnist()

#visualise some datas
par(mfrow = c(2,2))
for (i in 1:4){
  img <- t(apply(fashion$train$x[i,,], 2, rev))
  image(img,col=gray(0:255/255))
  print(fashion$train$y[i])
}
dev.off()

#Apprenstissage et validation des données
#Pretreatment datas
#Pretreament data
x_train <- fashion_mnist$train$x
y_train <- fashion_mnist$train$y
x_test <- fashion_mnist$test$x
y_test <- fashion_mnist$test$y
# reshape
x_train <- array_reshape(x_train, c(nrow(x_train), 784))
x_test <- array_reshape(x_test, c(nrow(x_test), 784))
# rescale
x_train <- x_train / 255
x_test <- x_test / 255
# output generation
y_train <- to_categorical(y_train, 10)
y_test <- to_categorical(y_test, 10)

model = keras_model_sequential()
#Training
model = layer_dense(model, units = 256, activation = "relu", input_shape = 784)
model = layer_dropout(model , rate = 0,4)
model = layer_dense(model, units = 128, activation = "sigmoid")
model = layer_dropout(model, rate = 0.1)
model = layer_dense(model, units = 10, activation = "softmax")

model=compile(model, loss = "categorical_crossentropy", optimizer = optimizer_rmsprop(), metrics = c("accuracy"))
history=fit(model, x_train, y_train, epochs = 20, batch_size = 128, validation_split = 0.2)

#visualisation
plot(history$metrics$acc)
points(history$metrics$val_acc,col="blue")

#Validation
y_pred=predict_classes(model, x_test)

#Evaluation
evaluate(model, x_test, y_test)
Table=table(y_pred, mnist$test$y)
sum(diag(Table))/10000



#******************************
#Exercice 3 : Critiques cinema
#******************************

imdb <- dataset_imdb(num_words = 500, maxlen = 100)

train_data=pad_sequences(imdb$train$x[1:4000], maxlen=100, value=0, padding="post")
test_data=pad_sequences(imdb$train$x[4001:5736], maxlen=100, value=0, padding="post")
train_labels=imdb$train$y[1:4000]
test_labels=imdb$train$y[4001:5736]

model <- keras_model_sequential()

model = layer_embedding(model, input_dim = 500, output_dim = 32, input_length=100)
model = layer_dropout(model, 0.25)
model = layer_lstm(model, units = 32)
model = layer_dense(model, units = 20, activation = "relu")
model = layer_dense(model, units = 1, activation = "softmax") # sortie

model = compile(model, optimizer = "rmsprop", loss = "binary_crossentropy", metrics = "accuracy")
result = fit(model, train_data, train_labels, epochs = 10, batch_size = 32, validation_split=0.1) # ou autre









