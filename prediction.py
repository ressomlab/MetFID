import numpy as np
import tensorflow as tf


def predict_fingerprint(binned_vector):
    """
    Given a binned vector, returns the predicted fingerprint by MetFID CNN
    model.
    :param binned_vector: a binned vector
    :return: predicted fingerprint
    """
    model = tf.keras.models.load_model('MetFID_CNN_40088_5618.h5')

    # for 1d CNN
    predicted = model.predict(np.array([binned_vector]), verbose=0)[0]

    # for 2d CNN
    # binned_vector = binned_vector.reshape(1, 35, 34, 1)
    # predicted = model.predict([binned_vector], verbose=0)[0]
    
    return np.array(list(map(lambda x: round(x), predicted)))
