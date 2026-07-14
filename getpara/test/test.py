import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import graphviz 
from sklearn.model_selection import train_test_split
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import LinearSVC
from sklearn.tree import DecisionTreeClassifier,export_graphviz

os.chdir('/Users/matsumi/teion/data/20220902/CH1_150mK_710uA_trig0.1_gain10_10kHz')
df = pd.read_csv('CH1/output/output_2.csv',index_col=0)
df = df.iloc[0:1000,:]

x = df.iloc[0:1000,2:7].values
y = df.iloc[0:1000,0].values
features = df.columns[2:7]
print(features)
x_train,x_test,y_train,y_test=train_test_split(x,y,test_size=0.3,random_state=0)

tree = DecisionTreeClassifier(max_depth=5,random_state=0)
tree.fit(x_train,y_train)
print("Training set accuracy: {:.5f}".format(tree.score(x_train,y_train)))
print("Test set accuracy: {:.5f}".format(tree.score(x_test,y_test)))

export_graphviz(tree,out_file="tree.dot",class_names=["1",'2'],feature_names=features,impurity=False,filled=True)

"""
clf = LinearSVC(C=1).fit(x_train,y_train)
#print("Training set prediction: {}".format(clf.predict(x_train)))
print("Training set accuracy: {:.5f}".format(clf.score(x_train,y_train)))
print("Test set accuracy: {:.5f}".format(clf.score(x_test,y_test)))
"""
"""
plt.plot(clf.coef_.T,'o')
plt.xticks(range(len(features)),features)
plt.xlabel('Feature')
plt.show()
"""

"""
neighbors_settings = range(1,11)
training_accuracy =[]
test_accuracy =[]

for n in neighbors_settings:
    clf = KNeighborsClassifier(n_neighbors=n)
    clf.fit(x_train,y_train)
    training_accuracy.append(clf.score(x_train,y_train))
    test_accuracy.append(clf.score(x_test,y_test))

print(len(clf.predict(x_train)))
print("Training set predictions: {}".format(clf.predict(x_train)))
#print("Training set accuracy: {}".format(training_accuracy))
#print("Test set accuracy: {}".format(test_accuracy))
"""
"""
plt.plot(neighbors_settings,training_accuracy,label="training accuracy")
plt.plot(neighbors_settings,test_accuracy,label="test accuracy")
plt.xlabel("n_neigbors")
plt.ylabel("Accuracy")
plt.legend()
plt.show()
"""

"""
pca =PCA(n_component=2)
svm = SVC()
sc = StandardScaler()
sc.fit(x_train)
x_scaler = sc.transform(x_train)
pca.fit(x_scaler)
x_pca = pca.transform(x_scaler)
svm.fit(x_pca,y_train)

"""

