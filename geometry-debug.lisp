(in-package #:org.shirakumo.fraf.convex-covering)

(defvar *material-number*)
(defun make-material (diffuse-factor &key (name (prog1
                                                    (format NIL "material~d" *material-number*)
                                                  (incf *material-number*)))
                                          alpha)
  (apply #'make-instance 'org.shirakumo.fraf.wavefront:material
         :name name
         :diffuse-factor diffuse-factor
         (when (and alpha (/= alpha 1))
           (list :transmission-factor (- 1 alpha)))))

(defvar *mesh-number*)
(defun make-triangle-mesh (vertices faces
                           &key (name (prog1
                                          (format NIL "mesh~d" *mesh-number*)
                                        (incf *mesh-number*)))
                                diffuse-factor
                                (alpha 1)
                                (material (make-material diffuse-factor :alpha alpha)))
  (make-instance 'org.shirakumo.fraf.wavefront:mesh :name name
                                                    :attributes '(:position)
                                                    :vertex-data vertices
                                                    :index-data faces
                                                    :material material))
