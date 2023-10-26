(in-package #:org.shirakumo.fraf.convex-covering)

;;; Debugging

(defvar *annotations*)
(defvar *annotation-number*)

(defun debug-line** (from direction &key (sample-count 32)
                                         (diffuse-factor #(.8 .8 .8)))
  (loop :with d = (v/ direction (float sample-count 1.0d0))
        :repeat sample-count
        :for v = from :then (v+ v d)
        :do (push (debug-cube v
                              (format nil "annotation~D" (incf *annotation-number*))
                              (make-material diffuse-factor)
                              :offset .001)
                  *annotations*)))

;;; Debug colors

(defvar *default-colors* #(#(1.0 0.0 0.0)
                           #(0.0 1.0 0.0)
                           #(0.0 0.0 1.0)
                           #(0.5 0.5 0.0)
                           #(0.0 0.5 0.5)
                           #(0.5 0.0 0.5)
                           #(1.0 1.0 1.0)
                           #(0.1 0.1 0.1)))

(defun make-color-generator (&key (colors *default-colors*))
  (let ((count (length colors)))
    (lambda (i)
      (aref colors (mod i count)))))

(defun color<-faces (faces)
  (flet ((face-value (index)
           (let ((a (aref faces (+ (* 3 index) 0)))
                 (b (aref faces (+ (* 3 index) 1)))
                 (c (aref faces (+ (* 3 index) 2))))
             (sxhash (logior (ash a 24) (ash b 12) c)))))
    (let ((hash (reduce #'logxor (alexandria:iota (/ (min 9 (length faces)) 3))
                        :key #'face-value)))
      (aref *default-colors* (mod hash (length *default-colors*))))))

;;; Debug geometry

(defun debug-cube-geometry (center offset)
  (let ((x (vx center))
        (y (vy center))
        (z (vz center)))
    (values (vector (+ x (- offset)) (+ y (- offset)) (+ z (- offset))
                    (+ x (- offset)) (+ y (- offset)) (+ z    offset )
                    (+ x (- offset)) (+ y    offset ) (+ z (- offset))
                    (+ x (- offset)) (+ y    offset ) (+ z    offset )
                    (+ x    offset ) (+ y (- offset)) (+ z (- offset))
                    (+ x    offset ) (+ y (- offset)) (+ z    offset )
                    (+ x    offset ) (+ y    offset ) (+ z (- offset))
                    (+ x    offset ) (+ y    offset ) (+ z    offset ))
            (vector 0 1 2  2 1 3
                    0 1 4  4 5 1
                    0 2 6  6 0 4
                    6 4 7  7 4 5
                    6 7 2  7 2 3
                    7 3 5  3 5 1))))

(defun debug-cube (center name material &key (offset .005))
  (multiple-value-bind (vertices faces) (debug-cube-geometry center offset)
    (make-instance 'org.shirakumo.fraf.wavefront:mesh
                   :name        name
                   :vertex-data vertices
                   :index-data  faces
                   :attributes  '(:position)
                   :material    material)))

(defun new-debug-line (from direction name material &key (offset .005)
                                                         (sample-count 32))
  (loop :with all-vertices = #()
        :with all-faces = #()
        :with d = (v/ direction (float sample-count 1.0d0))
        :repeat sample-count
        :for v = from :then (v+ v d)
        :for (vertices faces) = (multiple-value-list (debug-cube-geometry v offset))
        :do (setf all-faces    (concatenate 'vector all-faces (map 'vector (lambda (i)
                                                                             (+ (/ (length all-vertices) 3) i))
                                                                   faces))
                  all-vertices (concatenate 'vector all-vertices vertices))
        :finally (return (make-instance 'org.shirakumo.fraf.wavefront:mesh
                                        :name        name
                                        :vertex-data all-vertices
                                        :index-data  all-faces
                                        :attributes  '(:position)
                                        :material    material))))

(defun new-debug-line* (from to &rest args ;  &key offset sample-count
                                      )
  ; (declare (ignore offset sample-count))
  (apply #'new-debug-line from (v- to from) args))

(defun debug-line (from direction kind hull &key (sample-count (case kind
                                                                 (:boundary-edge 16)
                                                                 (t 32))))
  (loop :with d = (v/ direction (float sample-count 1.0d0))
        :repeat sample-count
        :for v = from :then (v+ v d)
        :do (push (cons v kind) (hull-annotations hull))))

(defun debug-line* (from to kind hull &rest args &key sample-count)
  (declare (ignore sample-count))
  (apply #'debug-line from (v- to from) kind hull args))

(defun debug-face (a b c kind hull &key (sample-count (case kind
                                                        (:boundary-edge 16)
                                                        (t 32))))
  (debug-line* a b kind hull :sample-count sample-count)
  (debug-line* b c kind hull :sample-count sample-count)
  (debug-line* c a kind hull :sample-count sample-count))

(defun debug-face* (vertices faces face-index kind hull &key (sample-count (case kind
                                                                             (:boundary-edge 16)
                                                                             (t 32))))
  (let ((a (manifolds:v vertices (+ (aref faces (+ (* 3 face-index) 0)))))
        (b (manifolds:v vertices (+ (aref faces (+ (* 3 face-index) 1)))))
        (c (manifolds:v vertices (+ (aref faces (+ (* 3 face-index) 2))))))
    (debug-face a b c kind hull :sample-count sample-count)))

;;; Utilities

(defun debug-filename (prefix i type)
  (format nil "/tmp/~A-~3,'0D.~A" prefix i type))

;;; Graph

(defclass decomposition ()
  ((%colors         :initarg  :colors
                    :reader   colors)
   (%order          :initarg  :order
                    :reader   order)
   (%bidirectinoal? :initarg  :both-directions?
                    :reader   both-directions?
                    :initform t)))

(defmethod cl-dot:graph-object-node ((graph decomposition) (object patch))
  (let* ((index (position object (order graph)))
         (color (when index
                  (color<-faces (hull-global-faces (patch-hull object)))
                  ; (funcall (colors graph) index)
                  ))
         (label (patch-debug-name object)))
    (make-instance 'cl-dot:node :attributes `(:label ,label
                                              :style :filled
                                              ,@(when color
                                                  `(:fillcolor ,(format nil "#~2,'0X~2,'0X~2,'0X"
                                                                        (+ 40 (floor (aref color 0) 1/215))
                                                                        (+ 40 (floor (aref color 1) 1/215))
                                                                        (+ 40 (floor (aref color 2) 1/215)))))))))

(defmethod cl-dot:graph-object-points-to ((graph decomposition) (object patch))
  (flet ((make-edge (link)
           (when (or (both-directions? graph)
                     (eq (patch-link-a link) object))
             (let* ((result     (patch-link-merge-result link))
                    (color      (cond ((eq link *winner*)
                                       "orange")
                                      ((typep result 'patch)
                                       "green")
                                      (t
                                       "red")))
                    (label      (when (typep result 'convex-hull)
                                  (format nil "~(~A~)" (hull-problem result))))
                    (width      (if (eq link *winner*) 5 nil))
                    (attributes `(:color ,color
                                  ,@(when width
                                      `(:penwidth ,width))
                                  ,@(when label
                                      `(:label ,label))
                                  ,@(when (not (both-directions? graph))
                                      '(:arrowhead :none)))))
               (make-instance 'cl-dot:attributed :object     (link-other-patch link object)
                                                 :attributes attributes)))))
    (remove nil (map 'list #'make-edge (patch-links object)))))

(defun graph-step (patches i &key (output-file (debug-filename "graph" i "png"))
                                  (colors      (make-color-generator)))
  (let* ((client (make-instance 'decomposition :colors colors :order (remove nil patches :key #'patch-hull)))
         (graph  (cl-dot:generate-graph-from-roots client patches)))
    (cl-dot:dot-graph graph output-file :format :png)))

;;; Render

(defun render-wavefront (input-file output-file &key (camera-position (vec 25 24 25)))
  (let ((camera-position (format nil "~A,~A,~A"
                                 (vx camera-position)
                                 (vy camera-position)
                                 (vz camera-position))))
    (inferior-shell:run `("f3d" "--output" ,output-file
                                "--camera-position" ,camera-position
                                ,input-file))))

(defun render-step (patches i &key (output-file (debug-filename "hulls" i "png"))
                                   colors
                                   highlight
                                   (camera-position (vec 25 24 25))
                                   (alpha .5))
  (let ((hulls           (remove nil (map 'vector #'patch-hull patches)))
        (object-filename (debug-filename "hulls" i "obj")))
    (export-hulls hulls object-filename :colors colors :alpha alpha :highlight highlight)
    (render-wavefront object-filename output-file :camera-position camera-position)))

;;; Visualization

(defun visualize-step (patches i &key highlight final)
  (when *debug-visualizations*
    (let ((patches (etypecase patches
                     (hash-table (alexandria:hash-table-keys patches))
                     (sequence   patches))))
      (let ((graph-filename (debug-filename "graph" i "png"))
            (image-filename (debug-filename "hulls" i "png"))
            (step-filename  (debug-filename "step" i "png")))
        (when (< (length patches) 64)
          (graph-step patches i :output-file graph-filename))
        (render-step patches i :output-file image-filename
                               :highlight   highlight
                               :colors      (when final
                                              (make-color-generator))
                               :alpha       nil #+no (unless final
                                                       .7))
        (when (< (length patches) 64)
          (inferior-shell:run `("montage" ,graph-filename ,image-filename
                                          "-tile" "2x1" "-geometry" "+0+0"
                                          ,step-filename)))))))

(defun visualize-problem (link i j)
  (return-from visualize-problem nil)
  (let* ((problem        (hull-problem (patch-link-merge-result link)))
         (debug-name1    (patch-debug-name (patch-link-a link) :format "~{~D-~D-~D~^_~}~@[__~]"))
         (debug-name2    (patch-debug-name (patch-link-b link) :format "~{~D-~D-~D~^_~}~@[__~]"))
         (base-name      (format nil "problem-~3,'0D/~A---~A-~A" i debug-name1 debug-name2 problem))
         (image-filename (debug-filename base-name j "png"))
                                        ; (step-filename  (debug-filename "step" i "png"))
         )
    (let ((hulls           (vector (patch-link-merge-result link)))
          (object-filename (debug-filename base-name j "obj")))
      (ensure-directories-exist image-filename)
      (export-hulls hulls object-filename :highlight (aref hulls 0))
      (inferior-shell:run `("f3d" "--output" ,image-filename
                                  "--camera-position" "5,4,5"
                                  ,object-filename)))))
