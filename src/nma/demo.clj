(ns nma.demo
  (:require [clojure.java.io :as jio]
            [incanter.core :as incanter]
            [nma.simple :as nma])
  (:import [javax.imageio ImageIO]
           [java.awt.image BufferedImage]))

(defn read-image-matrix [imagefile]
  (let [imagefile (jio/as-file imagefile)
        raster (.. (. ImageIO read imagefile) getData)
        left (. raster getMinX)
        width (. raster getWidth)
        right (+ left width)
        top (. raster getMinY)
        height (. raster getHeight)
        bottom (+ top height)
        xrange (range left right)
        yrange (range top bottom)]
    (incanter/matrix
     (for [row yrange]
       (for [column xrange]
         (let [pixel (. raster getPixel column row nil)]
           (/ (reduce + 0 (seq pixel)) (alength pixel))))))))

(defn write-image-matrix [imagefile M]
  (let [imagefile (jio/as-file imagefile)
        width (incanter/ncol M)
        height (incanter/nrow M)
        image (BufferedImage. width height BufferedImage/TYPE_BYTE_GRAY)
        raster (. image getRaster)]
    (doseq [row (range height) column (range width)]
      (. raster setPixel column row (float-array 1 (incanter/sel M row column))))
    (. ImageIO write image "png" imagefile)))


(defn image-demo
  "Will extract features from in-image-file, and write
   a reconstructed image to out-image-file. A bit useless,
   but a nice way to visualize the quality of a factorization,
   get a sense for how artifacts show up, etc."
  [in-image-file out-image-file features]
  (let [demo-src-image (read-image-matrix in-image-file)
        [W H] (nma/non-negative-factor-approx demo-src-image features
                                              :max-tries 100
                                              :chatty true)
        matrix-size (fn [M] (* (incanter/ncol M) (incanter/nrow M)))
        old-size (matrix-size demo-src-image)
        new-size (+ (matrix-size W) (matrix-size H))
        demo-reconstructed (incanter/mmult W H)]
    (println "Compressed image from " old-size " to " new-size)
    (write-image-matrix out-image-file demo-reconstructed)))
