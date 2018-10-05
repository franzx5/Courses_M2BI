library(keras)
#
# régression simple
#
# génération des données à partir modèle y=cos(x)+0.1*x^2
x=seq(-5,5,len=50)
ytrue=cos(x)+0.1*x*x
yobs=ytrue+rnorm(length(x),sd=0.2)
plot(x,ytrue,ty='l',col="red")
points(x,yobs, col="blue", pch="+")
# 
model = keras_model_sequential()
layer_dense(model, units=20, activation="relu", input_shape=1)
layer_dense(model, units=1, activation="linear", input_shape=4)
compile(model,loss = "mean_squared_error", optimizer = optimizer_rmsprop())
fit(model, x, yobs, epochs = 1500, batch_size = 10)
evaluate(model, x, yobs)
ypred=predict(model,x)
points(x,ypred, col="red")

#
# test classification simple: les iris
#
X=as.matrix(iris[,1:4])
y=as.integer(iris[,5])-1
Y=to_categorical(y)

model = keras_model_sequential()
layer_dense(model, units=5, activation="sigmoid", input_shape=4)
layer_dense(model, units=3, activation="softmax")
compile(model, loss = "categorical_crossentropy", optimizer = optimizer_rmsprop())
fit(model, X, Y, epochs = 1000, batch_size = 30)
evaluate(model, X, Y)
classes=predict_classes(model, X)
table(classes, y)

#
# test plus complet sur les iris
#
train=sample(1:150,100)
x_train = as.matrix(iris[train, 1:4])
y_train = iris[train, 5]
x_test = as.matrix(iris[-train, 1:4])
y_test = iris[-train,5]
y_train = to_categorical(as.integer(y_train)-1)
y_test = to_categorical(as.integer(y_test)-1)


model = keras_model_sequential() 
layer_dense(model, units = 5, activation = 'relu', input_shape = 4) 
layer_dense(model,units = 3, activation = 'softmax')
model
compile(model,loss = 'categorical_crossentropy',optimizer = optimizer_rmsprop(),
  metrics = 'accuracy')
fit(model, x_train, y_train, epochs = 500, batch_size = 50, validation_split = 0.2)
evaluate(model, x_test, y_test)
classes=predict_classes(model, x_train)
table(classes, iris[train,5])
classes=predict_classes(model, x_test)
table(classes, iris[-train,5])

#
# données MNIST 60000 images 28x28 et labels
#
mnist <- dataset_mnist()
k=5
img <- t(apply(mnist$train$x[5,,], 2, rev))
image(img,col=gray(0:255/255))
mnist$train$y[k]

#
mnist <- dataset_mnist()
x_train <- mnist$train$x
y_train <- mnist$train$y
x_test <- mnist$test$x
y_test <- mnist$test$y
# reshape
x_train <- array_reshape(x_train, c(nrow(x_train), 784))
x_test <- array_reshape(x_test, c(nrow(x_test), 784))
# rescale
x_train <- x_train / 255
x_test <- x_test / 255
# output generation
y_train <- to_categorical(y_train, 10)
y_test <- to_categorical(y_test, 10)


model=keras_model_sequential() 
#input_shape obligatoire comme premiere couche, dim. output=256
model = layer_dense(model,units = 256, activation = "relu", input_shape = 784)
#input_shape 
model = layer_dropout(model,rate = 0.4)
model = layer_dense(model,units = 128, activation = "relu") 
model = layer_dropout(model,rate = 0.3) 
model = layer_dense(model,units = 10, activation = "softmax")

model=compile(model, loss = "categorical_crossentropy", optimizer = optimizer_rmsprop(), metrics = c("accuracy"))

history=fit(model, x_train, y_train, epochs = 30, batch_size = 128, validation_split = 0.2)

plot(history$metrics$acc)
points(history$metrics$val_acc,col="blue")

y_pred=predict_classes(model, x_test)

evaluate(model, x_test, y_test)
Table=table(y_pred, mnist$test$y)
sum(diag(Table))/10000

#
# données textuelles
#
imdb <- dataset_imdb(num_words = 500, maxlen = 100)
train_data=pad_sequences(imdb$train$x[1:4000], maxlen=100, value=0, padding="post")
test_data=pad_sequences(imdb$train$x[4001:5736], maxlen=100, value=0, padding="post")
train_labels=imdb$train$y[1:4000]
test_labels=imdb$train$y[4001:5736]
model <- keras_model_sequential()
model = layer_embedding(model, input_dim = 500, output_dim = 32, input_length=100)
model = layer_dropout() 
# model rate = 0.25
# units 32, units 256
model = layer_lstm() # avec ou sans
model = layer_dense() # couche cachée
model = layer_dense() # sortie
model = compile(model, optimizer = "rmsprop", loss = "binary_crossentropy", metrics = "accuracy")
result = fit(model, train_data, train_labels, epochs = 10, batch_size = 32, validation_split=0.1) # ou autre
