(in-package #:org.shirakumo.fraf.convex-covering)

;;; Variables

(defvar *debug-output* nil)

(defvar *debug-visualizations* nil)

(defun debug-visualizations-p (i)
  (and *debug-visualizations*
       (zerop (mod i *debug-visualizations*))))

;;; Debug output

(defun d (format-control &rest format-arguments)
  (when *debug-output*
    (apply #'format *trace-output* format-control format-arguments)))

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

;;; Debug geometry

(defun debug-cube (center name material &key (offset .005))
  (let ((x (vx center))
        (y (vy center))
        (z (vz center)))
    (make-instance 'org.shirakumo.fraf.wavefront:mesh
                   :name name
                   :vertex-data (vector (+ x (- offset)) (+ y (- offset)) (+ z (- offset))
                                        (+ x (- offset)) (+ y (- offset)) (+ z    offset )
                                        (+ x (- offset)) (+ y    offset ) (+ z (- offset))
                                        (+ x (- offset)) (+ y    offset ) (+ z    offset )
                                        (+ x    offset ) (+ y (- offset)) (+ z (- offset))
                                        (+ x    offset ) (+ y (- offset)) (+ z    offset )
                                        (+ x    offset ) (+ y    offset ) (+ z (- offset))
                                        (+ x    offset ) (+ y    offset ) (+ z    offset ))
                   :index-data (vector 0 1 2  2 1 3
                                       0 1 4  4 5 1
                                       0 2 6  6 0 4
                                       6 4 7  7 4 5
                                       6 7 2  7 2 3
                                       7 3 5  3 5 1)
                   :attributes '(:position)
                   :material material)))

(defun debug-line (from direction kind hull &key (sample-count 32))
  (loop :with d = (v/ direction (float sample-count 1.0d0))
        :repeat sample-count
        :for v = from :then (v+ v d)
        :do (push (cons v kind) (annotations hull))))

(defun debug-line* (from to kind hull &rest args &key sample-count)
  (declare (ignore sample-count))
  (apply #'debug-line from (v- to from) kind hull args))

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
                  (funcall (colors graph) index)))
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
                    (label      (unless (typep result 'patch)
                                  (format nil "~(~A~)" (problem result))))
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

(defun render-step (patches i &key (output-file (debug-filename "hulls" i "png"))
                                   (colors      (make-color-generator))
                                   highlight)
  (let ((hulls           (remove nil (map 'vector #'patch-hull patches)))
        (object-filename (debug-filename "hulls" i "obj")))
    (org.shirakumo.fraf.convex-covering.test::export-hulls
     hulls object-filename :colors colors :highlight highlight)
    (inferior-shell:run `("f3d" "--output" ,output-file
                                "--camera-position" "5,4,5"
                                ,object-filename))))

;;; Visualization

(defun visualize-step (patches i &key highlight)
  (when *debug-visualizations*
    (let ((patches (etypecase patches
                     (hash-table (alexandria:hash-table-keys patches))
                     (sequence   patches))))
      (let ((graph-filename (debug-filename "graph" i "png"))
            (image-filename (debug-filename "hulls" i "png"))
            (step-filename  (debug-filename "step" i "png")))
        (when (< (length patches) 64)
          (graph-step patches i :output-file graph-filename))
        (render-step patches i :output-file image-filename :highlight highlight)
        (when (< (length patches) 64)
          (inferior-shell:run `("montage" ,graph-filename ,image-filename
                                          "-tile" "2x1" "-geometry" "+0+0"
                                          ,step-filename)))))))

(defun visualize-problem (link i j)
  (return-from visualize-problem nil)
  (let* ((problem        (problem (patch-link-merge-result link)))
         (debug-name1    (patch-debug-name (patch-link-a link) :format "~{~D-~D-~D~^_~}~@[__~]"))
         (debug-name2    (patch-debug-name (patch-link-b link) :format "~{~D-~D-~D~^_~}~@[__~]"))
         (base-name      (format nil "problem-~3,'0D/~A---~A-~A" i debug-name1 debug-name2 problem))
         (image-filename (debug-filename base-name j "png"))
                                        ; (step-filename  (debug-filename "step" i "png"))
         )
    (let ((hulls           (vector (patch-link-merge-result link)))
          (object-filename (debug-filename base-name j "obj")))
      (ensure-directories-exist image-filename)
      (org.shirakumo.fraf.convex-covering.test::export-hulls
       hulls object-filename :highlight (aref hulls 0))
      (inferior-shell:run `("f3d" "--output" ,image-filename
                                  "--camera-position" "5,4,5"
                                  ,object-filename)))))
