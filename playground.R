#the data set
movie_features = matrix(runif(100,1,10),nrow = 50,ncol = 2)
user_features = matrix(runif(20,1,5),nrow = 10, ncol = 2)
training_movies = movie_features[c(1:40),]
predict_movies = movie_features[c(41:50),]

training_ratings = training_movies %*% t(user_features)
actual_ratings = predict_movies %*% t(user_features)
all_ratings = movie_features %*% t(user_features)

movie_df = data.frame(feat = movie_features)
rownames(movie_df) = paste("movie",c(1:dim(movie_df)[1]),sep = "")

user_df = data.frame(feat = user_features)
rownames(user_df) = paste("user",c(1:dim(user_df)[1]),sep = "")

ratings_df = data.frame(users = all_ratings)
rownames(ratings_df) = paste("movie",c(1:dim(movie_df)[1]),sep = "")

incomp_df = ratings_df
incomp_df[floor(runif(10,1,50)),c(1:3)] = -1
incomp_df[floor(runif(10,1,50)),c(4:7)] = -1
incomp_df[floor(runif(10,1,50)),c(8:10)] = -1

learned_user_features = matrix(runif(20,1,2),nrow = 10, ncol = 2)
learned_movie_features = matrix(0,nrow= dim(movie_features)[1], ncol = dim(movie_features)[2])
#ridge regression to find movie feature vectors

for (i in c(1:20)){
  for (i in c(1:dim(ratings_df)[1])){
    ridge = lm.ridge(t(data.matrix(incomp_df[i,incomp_df[i,] != -1])) ~ feat.1 + feat.2 + 0, data=data.frame( feat = learned_user_features[incomp_df[i,] != -1,]), lambda=0)
    learned_movie_features[i,] = coef(ridge)
    #incomp_df[i,incomp_df[i,] == -1] = coef(ridge)[1] +  user_features[incomp_df[i,] == -1,]  %*% data.matrix(coef(ridge)[c(2:3)])
  }
  
  
  #ridge regression to find user feature vector
  for (i in c(1:dim(ratings_df)[2])){
    ridge = lm.ridge(incomp_df[incomp_df[,i] != -1,i] ~ feat.1 + feat.2 + 0 , data=data.frame( feat = learned_movie_features[incomp_df[,i] != -1,]), lambda=0)
    learned_user_features[i,] = coef(ridge)
    #incomp_df[incomp_df[,i] == -1,i] = coef(ridge)[1] +  movie_features[incomp_df[,i] == -1,]  %*% data.matrix(coef(ridge)[c(2:3)])
  }
  print(sum(abs(all_ratings - learned_movie_features %*% t(learned_user_features))))
}



#lasso regression
#come back to here

#simple regression
if (FALSE){
  model = lm(training_ratings[,1] ~ feat.1 + feat.2, data = data.frame(feat = training_movies))
  predictions = predict(model,data.frame(feat = predict_movies))
  sum(predictions - actual_ratings[,1])
}
