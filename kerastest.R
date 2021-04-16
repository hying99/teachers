
data <- matrix(runif(1000*100), nrow = 1000, ncol = 100)
labels <- matrix(round(runif(1000, min = 0, max = 1)), nrow = 1000, ncol = 1)

model <- keras_model_sequential()
K <- backend()
 metric_mean_pred <- function(y_true, y_pred) {
  K$dot(y_true,y_pred)
   }

model %>%
  layer_dense(units = 32, activation = 'relu', input_shape = c(100)) %>%
  layer_dense(units = 1, activation = 'sigmoid') %>%
  compile(    optimizer = 'rmsprop',
    loss = 'binary_crossentropy',
    metrics = c(metric_binary_crossentropy)  )

model %>% fit(data, labels, epochs=8, batch_size=32)




























# mnist <- dataset_mnist()  #导入mnist数据集
# x_train <- mnist$train$x   #训练集的自变量
# y_train <- mnist$train$y   #训练集的因变量
# x_test <- mnist$test$x
# y_test <- mnist$test$y
# #改变数据形状，矩阵转向量
# dim(x_train) <- c(nrow(x_train), 784)
# dim(x_test) <- c(nrow(x_test), 784)
# # 归一化
# x_train <- x_train / 255
# x_test <- x_test / 255
# #调整输出的形式（因变量）
# y_train <- to_categorical(y_train, 10)
# y_test <- to_categorical(y_test, 10)
# model_total=list()
# model <- keras_model_sequential()
# model %>%
#   layer_dense(units = 256, activation = 'relu', input_shape = c(784)) %>%
#   layer_dropout(rate = 0.4) %>%
#   layer_dense(units = 128, activation = 'relu') %>%
#   layer_dropout(rate = 0.3) %>%
#   layer_dense(units = 10, activation = 'softmax')
# #自定义metric函数
# K <- backend()
# metric_mean_pred <- function(y_true, y_pred) {
#   K$mean(y_pred) 
# }
# 
# 
# model %>% compile(
#   loss = loss_mean_squared_error,
#   optimizer = optimizer_adagrad(),
#   #metrics = c(metric_categorical_crossentropy)
#   metrics = c(metric_mean_pred)
# )
# model_total[[1]]=model
# 
# history <- model %>% fit(
#   x_train, y_train,
#   epochs = 5, batch_size = 256,
#   #validation_split = 0.2
#   validation_data = list(x_test,y_test))

# history.1 <- model_total[[1]] %>% fit(
#   x_train, y_train, 
#   epochs = 3, batch_size = 300, 
#   validation_split = 0.2)