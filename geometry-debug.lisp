(in-package #:org.shirakumo.fraf.convex-covering)

(defvar *material-number*)
(defun make-material (diffuse-factor &key (name (prog1
                                                    (format nil "material~D" *material-number*)
                                                  (incf *material-number*)))
                                          alpha)
  (apply #'make-instance 'org.shirakumo.fraf.wavefront:material
         :name name
         :diffuse-factor diffuse-factor
         (when (and alpha (/= alpha 1))
           (list :transmission-factor (- 1 alpha)))))

(defvar *mesh-number*)
(defun make-triangle-mesh (vertices faces &key (name (prog1
                                                         (format nil "mesh~D" *mesh-number*)
                                                       (incf *mesh-number*)))
                                               diffuse-factor
                                               (alpha 1)
                                               (material (make-material diffuse-factor :alpha alpha)))
  (make-instance 'org.shirakumo.fraf.wavefront:mesh
                 :name name
                 :attributes '(:position)
                 :vertex-data vertices
                 :index-data faces
                 :material material))

;;; Coplanar cases
(defun test-1 ()
  (let ((*material-number* 1)
        (*mesh-number* 1)
        (*annotations* '())
        (*annotation-number* 0)
        (*debug-output* t))
    (multiple-value-bind (vertices-u faces-u) (gen-triangle)
      (multiple-value-bind (vertices-v faces-v) (gen-triangle)
        (multiple-value-bind (vertices-w faces-w) (gen-triangle)
          (multiple-value-bind (vertices-x faces-x) (gen-triangle)
            (multiple-value-bind (vertices-y faces-y) (gen-triangle)
              (progn ; multiple-value-bind (vertices-z faces-z) (gen-triangle)
                ;; Coplanar intersecting
                (let ((c (v/ (reduce #'v+ vertices-u) 3)))
                  (map-into vertices-v (lambda (v) (v+ c (v* (v- v c) (+ .7 (random .6))))) vertices-u))
                ;; Coplanar non-intersecting
                (let ((d (v- (aref vertices-u 1) (aref vertices-u 0))))
                  (map-into vertices-w (lambda (v) (v+ v (v* d 2))) vertices-u))
                ;; Coplanar one-edge intersection other vertex outside
                (let ((d (v- (aref vertices-u 2) (aref vertices-u 0))))
                  (setf (aref vertices-x 0) (aref vertices-u 0)
                        (aref vertices-x 1) (aref vertices-u 1)
                        (aref vertices-x 2) (v+ (aref vertices-u 0) (v* d -2))))
                ;; Coplanar one-edge intersection other vertex inside
                (let ((c (v/ (v+ (aref vertices-u 0) (aref vertices-u 1) (aref vertices-u 2)) 3)))
                  (setf (aref vertices-y 0) (aref vertices-u 0)
                        (aref vertices-y 1) (aref vertices-u 1)
                        (aref vertices-y 2) c))

                (let ((intersect-uv-p (triangles-intersect-p (aref vertices-u 0)
                                                             (aref vertices-u 1)
                                                             (aref vertices-u 2)
                                                             (aref vertices-v 0)
                                                             (aref vertices-v 1)
                                                             (aref vertices-v 2)))
                      (intersect-uw-p (triangles-intersect-p (aref vertices-u 0)
                                                             (aref vertices-u 1)
                                                             (aref vertices-u 2)
                                                             (aref vertices-w 0)
                                                             (aref vertices-w 1)
                                                             (aref vertices-w 2)))
                      (intersect-ux-p (triangles-intersect-p (aref vertices-u 0)
                                                             (aref vertices-u 1)
                                                             (aref vertices-u 2)
                                                             (aref vertices-x 0)
                                                             (aref vertices-x 1)
                                                             (aref vertices-x 2)))
                      (intersect-uy-p (triangles-intersect-p (aref vertices-u 0)
                                                             (aref vertices-u 1)
                                                             (aref vertices-u 2)
                                                             (aref vertices-y 0)
                                                             (aref vertices-y 1)
                                                             (aref vertices-y 2)))
                      #+no (intersect-uz-p (triangles-intersect-p (aref vertices-u 0)
                                                                  (aref vertices-u 1)
                                                                  (aref vertices-u 2)
                                                                  (aref vertices-z 0)
                                                                  (aref vertices-z 1)
                                                                  (aref vertices-z 2))))
                  (let* ((object-file #P"/tmp/triangle-intersection.obj")
                         (output-file #P"/tmp/triangle-intersection.png")
                         (objects     (append (list (make-triangle-mesh (flatten-vertices vertices-u) faces-u
                                                                        :diffuse-factor #(.7 .5 .5))
                                                    #+no (make-triangle-mesh (flatten-vertices (offset-vertices-along-normal vertices-u bias)) faces-u
                                                                             :diffuse-factor #(.7 1 1)
                                                                             :alpha .5)
                                                    #+no (make-triangle-mesh (flatten-vertices (offset-vertices-along-normal vertices-u (- bias))) faces-u
                                                                             :diffuse-factor #(1 1 .7)
                                                                             :alpha .5)
                                                    (make-triangle-mesh (flatten-vertices vertices-v) faces-v
                                                                        :diffuse-factor (if intersect-uv-p
                                                                                            #(0 1 0)
                                                                                            #(.6 .7 .6)))
                                                    #+no (make-triangle-mesh (flatten-vertices vertices-w) faces-w
                                                                             :diffuse-factor (if intersect-uw-p
                                                                                                 #(0 0 1)
                                                                                                 #(.6 .6 .7)))
                                                    (make-triangle-mesh (flatten-vertices vertices-x) faces-x
                                                                        :diffuse-factor (if intersect-ux-p
                                                                                            #(1 0 1)
                                                                                            #(.7 .6 .7)))
                                                    (make-triangle-mesh (flatten-vertices vertices-y) faces-y
                                                                        :diffuse-factor (if intersect-uy-p
                                                                                            #(0 1 1)
                                                                                            #(.6 .7 .7)))
                                                    #+no (make-triangle-mesh (flatten-vertices vertices-z) faces-z
                                                                             :diffuse-factor (if intersect-uz-p
                                                                                                 #(1 1 0)
                                                                                                 #(.7 .7 .6))))
                                              *annotations*)))
                    (org.shirakumo.fraf.wavefront:serialize objects object-file :if-exists :supersede)
                    (render-wavefront object-file output-file :camera-position (dvec 2 2 2))))))))))))
