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

;;; Helper functions

(defun p~ (point1 point2 &key (eps 1d-6))
  (<= (vsqrdistance point1 point2) eps))

;;; Intersections via Separating Axis Theorem

(defun intervals-intersect-p (i1-min i1-max i2-min i2-max &key (threshold 1d-6))
  (declare (type double-float i1-min i1-max i2-min i2-max)
           (optimize (speed 3)))
  ;; case 1 | s1 e1 s2 e2
  ;; case 2 | s1 s2 e1 e2 *
  ;; case 3 | s1 s2 e2 e1 *
  ;; case 4 | s2 s1 e1 e2 *
  ;; case 5 | s2 s1 e2 e1 *
  ;; case 6 | s2 e2 s1 e1
  #+no (or (<= i1-min i2-min i1-max)
           (<= i2-min i1-min i2-max))
  (cond ((<= (abs (- (min i1-min i2-min) (max i1-max i2-max))) threshold)
         (values t :identical))
        ((<= i1-min i2-min i1-max i2-max) ; case 2
         (d "~4@TOverlap ~A~%"
            (abs (- i2-min i1-max)))
         (values t (if (< (abs (- i2-min i1-max)) threshold) :touching :penetrating)))
        ((<= i1-min i2-min i2-max i1-max) ; case 3
         (d "~4@TI₂ ⊂ I₁ (case 3) Overlaps ~A ~A~%"
            (abs (- i1-min i2-min)) ; cases 3 4
            (abs (- i1-max i2-max)) ; cases 3 4
            )
         (values t :penetrating
                 #+no (if (< (min (abs (- i1-min i2-min)) (abs (- i1-max i2-max))) threshold)
                       :touching
                       :penetrating)))
        ((<= i2-min i1-min i1-max i2-max) ; case 4
         (d "~4@TI₁ ⊂ I₂ (case 4) Overlaps ~A ~A~%"
            (abs (- i1-min i2-min)) ; cases 3 4
            (abs (- i1-max i2-max)) ; cases 3 4
            )
         (values t :penetrating
                 #+no (if (< (min (abs (- i1-min i2-min)) (abs (- i1-max i2-max))) threshold) :touching :penetrating)))
        ((<= i2-min i1-min i2-max i1-max) ; case 5
         (d "~4@TOverlap ~A~%"
            (abs (- i1-min i2-max)))
         (values t (if (< (abs (- i1-min i2-max)) threshold) :touching :penetrating)))
        #+no ((not (or (> (- i2-min i1-max) threshold)
                  (> (- i1-min i2-max) threshold)))
         (d "~4@TOverlaps ~A ~A ~A ~A~%"
            (abs (- i2-min i1-max)) ; case 2
            (abs (- i1-min i2-min)) ; cases 3 4
            (abs (- i1-max i2-max)) ; cases 3 4
            (abs (- i1-min i2-max))) ; case 5
         t)
        (t ; no intersection
         nil)))

(defun triangles-intersect-old-p (u1 u2 u3 v1 v2 v3 &key (bias .0001))
  (flet ((separating-axis-p (axis bias &optional note)
           (let* ((du1    (+ (v. axis u1) bias))
                  (du2    (+ (v. axis u2) bias))
                  (du3    (+ (v. axis u3) bias))
                  (dv1    (v. axis v1))
                  (dv2    (v. axis v2))
                  (dv3    (v. axis v3))
                  (iu-min (min du1 du2 du3)) ; interval u min
                  (iu-max (max du1 du2 du3)) ; interval u max
                  (iv-min (min dv1 dv2 dv3)) ; interval v min
                  (iv-max (max dv1 dv2 dv3)))
             (d "Axis~@[ ~A~] ~A~&~
                 Bias ~A~&~
                 ~2@T=> {~5,2F ~5,2F ~5,2F} {~5,2F ~5,2F ~5,2F}~&~
                 ~2@T=> [~5,2F, ~5,2F] [~5,2F, ~5,2F]~%"
                note axis bias
                du1 du2 du3 dv1 dv2 dv3
                iu-min iu-max iv-min iv-max)
             (let ((result (intervals-intersect-p iu-min iu-max iv-min iv-max)))
               (format *trace-output* "~2@T~A~%" (if result "?" "separating"))
               (not result)))))
    (let* ((normal-u (vunit (vc (v- u2 u1) (v- u3 u1))))
           (n12      (vunit (vc normal-u (v- u2 u1)))) ; TODO compute when needed
           (n23      (vunit (vc normal-u (v- u3 u2))))
           (n31      (vunit (vc normal-u (v- u1 u3)))))

      (apply #'debug-line** (v/ (v+ u1 u2 u3) 3) (v* normal-u .1)
             (when (or (separating-axis-p normal-u      bias)
                       (separating-axis-p (v- normal-u) (- bias)))
               (list :diffuse-factor #(1 0 0))))
      (apply #'debug-line** (v/ (v+ u1 u2) 2)    (v* n12 .05)
             (when (separating-axis-p n12 0)
               (list :diffuse-factor #(1 0 0))))
      (apply #'debug-line** (v/ (v+ u2 u3) 2)    (v* n23 .05)
             (when (separating-axis-p n23 0)
               (list :diffuse-factor #(1 0 0))))
      (apply #'debug-line** (v/ (v+ u3 u1) 2)    (v* n31 .05)
             (when (separating-axis-p n31 0)
               (list :diffuse-factor #(1 0 0))))

      (and (not (separating-axis-p normal-u      bias     "normal"))
           (not (separating-axis-p (v- normal-u) (- bias) "-normal"))
           (not (separating-axis-p n12           0))
           (not (separating-axis-p n23           0))
           (not (separating-axis-p n31           0))))))

(defun triangles-intersect-p (u1 u2 u3 v1 v2 v3 &key (bias .0001) (threshold 1d-6))
  (d "U~2@T~A~%~2@T~A~%~2@T~A~%V~2@T~A~%~2@T~A~%~2@T~A~%" u1 u2 u3 v1 v2 v3)
  (flet ((separating-axis-p (axis bias &optional note)
           (let* ((du1    (+ (v. axis u1) bias))
                  (du2    (+ (v. axis u2) bias))
                  (du3    (+ (v. axis u3) bias))
                  (dv1    (v. axis v1))
                  (dv2    (v. axis v2))
                  (dv3    (v. axis v3))
                  (iu-min (min du1 du2 du3)) ; interval u min
                  (iu-max (max du1 du2 du3)) ; interval u max
                  (iv-min (min dv1 dv2 dv3)) ; interval v min
                  (iv-max (max dv1 dv2 dv3)))
             (d "~2@TAxis~@[ ~A~] ~A~&~
                 ~4@T=> {~5,2F ~5,2F ~5,2F} {~5,2F ~5,2F ~5,2F}~&~
                 ~4@T=> [~5,2F, ~5,2F] [~5,2F, ~5,2F]~%"
                note axis
                du1 du2 du3 dv1 dv2 dv3
                iu-min iu-max iv-min iv-max)
             (multiple-value-bind (intersectp class)
                 (intervals-intersect-p iu-min iu-max iv-min iv-max
                                        :threshold threshold)
               (d "~4@T=> ~A~@[ | ~A~]~%" (if intersectp "not separating" "separating") class)
               (values (not intersectp) class)))))
    (d "--------------------~%")
    (multiple-value-bind (result constellation-class contact-class)
        (let ((normal-u (vunit (vc (v- u2 u1) (v- u3 u1)))))
          (d "Normal~%")
          (multiple-value-bind (normal-separating-p normal-class)
              (separating-axis-p normal-u 0)
            (d "~2@T=> ~:[not separating~;separating~]~@[ | ~A~]~%"
               normal-separating-p normal-class)
            (cond (normal-separating-p
                   (d "normal separating~%")
                   nil)
                  ((eq normal-class :identical)
                   (d "Edge normals (coplanar case)~%")
                   (block nil
                     (let ((contact-class :penetrating))
                       (flet ((test-axis (e1 e2)
                                (declare (type dvec3 e1 e2))
                                (let ((edge-normal (vunit (vc normal-u (v- e2 e1)))))
                                  (multiple-value-bind (separatingp class)
                                      (separating-axis-p edge-normal 0)
                                    (cond (separatingp
                                           (return (values nil :coplanar)))
                                          (t
                                           (when (and (eq contact-class :penetrating)
                                                      (eq class :touching))
                                             (setf contact-class class))))))))
                         (test-axis u1 u2)
                         (test-axis u2 u3)
                         (test-axis u3 u1)
                         (test-axis v1 v2)
                         (test-axis v2 v3)
                         (test-axis v3 v1)
                         (values t :coplanar contact-class))))

                   #+old(let ((n12 (vunit (vc normal-u (v- u2 u1))))
                              (n23 (vunit (vc normal-u (v- u3 u2))))
                              (n31 (vunit (vc normal-u (v- u1 u3)))))

                          (cond ((and (not (separating-axis-p n12 0))
                                      (not (separating-axis-p n23 0))
                                      (not (separating-axis-p n31 0)))
                                 (values t :coplanar))
                                (t
                                 (values nil :coplanar)))))
                  (t
                   (d "General edge pairs~%")
                   (let ((contact-class normal-class))
                     (if (flet ((test-axis (u1 u2 v1 v2)
                                  (let ((c (vc (v- u2 u1) (v- v2 v1))))
                                    (unless (< (vsqrlength c) threshold)
                                      (let ((c* (vunit c)))
                                        (multiple-value-bind (separatingp class)
                                            (separating-axis-p c* 0)
                                          (when (and (not separatingp)
                                                     (eq contact-class :penetrating)
                                                     (eq class :touching))
                                            (setf contact-class class))
                                          (when (boundp '*annotation-number*)
                                            (apply #'debug-line** u1 (v* c* .1)
                                                   (when separatingp
                                                     (list :diffuse-factor #(1 0 0)))))
                                          separatingp))))))
                           (prog1
                               (loop :named outer
                                     :for (a1 a2) :on (list u1 u2 u3 u1)
                                     :while a2
                                     :do (loop :for (b1 b2) :on (list v1 v2 v3 v1)
                                               :while b2
                                               :when (test-axis a1 a2 b1 b2)
                                               :do (return-from outer t)))))
                         nil
                         (values t nil contact-class)))))))
      (d "=> ~A~:*~:[~; ~A ~A~]~%"
         result
         (or constellation-class :general)
         contact-class)
      (d "--------------------~%")
      (values result constellation-class contact-class))))

(defun line-intersects-triangle-p (u1 u2 u3 l1 l2 &key (eps 1d-6))
  )

(defun point-on-triangle-p (u1 u2 u3 v &key (threshold .00001))
  (flet ((separating-axis-p (axis bias &optional note)
           (let* ((du1    (+ (v. axis u1) bias))
                  (du2    (+ (v. axis u2) bias))
                  (du3    (+ (v. axis u3) bias))
                  (dv     (v. axis v))
                  (iu-min (min du1 du2 du3)) ; interval u min
                  (iu-max (max du1 du2 du3)) ; interval u max
                  )
             (d "Axis~@[ ~A~] ~A~&~
                 ~2@T=> {~5,2F ~5,2F ~5,2F} {~5,2F ~5,2F ~5,2F}~&~
                 ~2@T=> [~5,2F, ~5,2F] [~5,2F, ~5,2F]~%"
                note axis bias
                du1 du2 du3 dv dv dv
                iu-min iu-max dv dv)
             (let ((result (intervals-intersect-p iu-min iu-max dv dv
                                                  :threshold threshold)))
               (d "~2@T~A~%" (if result "not separating" "separating"))
               (not result)))))
    (let* ((normal-u (vunit (vc (v- u2 u1) (v- u3 u1))))
           (n12      (vunit (vc normal-u (v- u2 u1)))) ; TODO compute when needed
           (n23      (vunit (vc normal-u (v- u3 u2))))
           (n31      (vunit (vc normal-u (v- u1 u3)))))
      (d "--------------------~%")
      (prog1
          (and (< (abs (v. normal-u (v- v u1))) threshold)
               (not (separating-axis-p n12 0))
               (not (separating-axis-p n23 0))
               (not (separating-axis-p n31 0)))
        (d "--------------------~%")))))

;;; Test

(defun gen-triangle (&key (gen-coord  (lambda ()
                                        (- (random 1.0d0) 0.5d0)))
                          (gen-vertex (lambda ( ; gen-coord
                                               )
                                        (flet ((coord ()
                                                 (funcall gen-coord)))
                                          (dvec (coord) (coord) (coord))))))
  (let ((a (funcall gen-vertex))
        (b (funcall gen-vertex))
        (c (funcall gen-vertex)))
    (values (vector a b c) (vector 0 1 2))))

(defun flatten-vertices (vertices)
  (let ((result (make-array (* 3 (length vertices))))
        (i      0))
    (map nil (lambda (vertex)
               (setf (aref result (+ i 0)) (vx vertex)
                     (aref result (+ i 1)) (vy vertex)
                     (aref result (+ i 2)) (vz vertex))
               (incf i 3))
         vertices)
    result))

(defun offset-vertices (vertices direction))

(defun offset-vertices-along-normal (vertices offset)
  (let* ((a      (aref vertices 0))
         (b      (aref vertices 1))
         (c      (aref vertices 2))
         (normal (vc (v- b a) (v- c a))))
    (vector (v+ a (v* normal offset))
            (v+ b (v* normal offset))
            (v+ c (v* normal offset)))))

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

(defun test-2 ()
  (let ((*material-number* 1)
        (*mesh-number* 1)
        (*annotations* '())
        (*annotation-number* 0)
        (*debug-output* t)

        (bias .125))
    (multiple-value-bind (vertices-u faces-u) (gen-triangle)
      (multiple-value-bind (vertices-v faces-v) (gen-triangle)
        (multiple-value-bind (vertices-w faces-w) (gen-triangle)
          (multiple-value-bind (vertices-x faces-x) (gen-triangle)
            (multiple-value-bind (vertices-y faces-y) (gen-triangle)
              (multiple-value-bind (vertices-z faces-z) (gen-triangle)
                (let ((c (v/ (reduce #'v+ vertices-u) 3)))
                  (map-into vertices-v (lambda (v) (v+ c (v* (v- v c) (+ .7 (random .6))))) vertices-u))
                (let ((d (v- (aref vertices-u 1) (aref vertices-u 0))))
                  (map-into vertices-w (lambda (v) (v+ v (v* d 2))) vertices-u))
                (setf (aref vertices-x 0) (aref vertices-u 0)
                      (aref vertices-x 1) (aref vertices-u 1))
                (setf (aref vertices-y 0) (aref vertices-u 0))
                (let ((intersect-uv-p (triangles-intersect-p (aref vertices-u 0)
                                                             (aref vertices-u 1)
                                                             (aref vertices-u 2)
                                                             (aref vertices-v 0)
                                                             (aref vertices-v 1)
                                                             (aref vertices-v 2)
                                                             :bias bias))
                      (intersect-uw-p (triangles-intersect-p (aref vertices-u 0)
                                                             (aref vertices-u 1)
                                                             (aref vertices-u 2)
                                                             (aref vertices-w 0)
                                                             (aref vertices-w 1)
                                                             (aref vertices-w 2)
                                                             :bias bias))
                      (intersect-ux-p (triangles-intersect-p (aref vertices-u 0)
                                                             (aref vertices-u 1)
                                                             (aref vertices-u 2)
                                                             (aref vertices-x 0)
                                                             (aref vertices-x 1)
                                                             (aref vertices-x 2)
                                                             :bias bias))
                      (intersect-uy-p (triangles-intersect-p (aref vertices-u 0)
                                                             (aref vertices-u 1)
                                                             (aref vertices-u 2)
                                                             (aref vertices-y 0)
                                                             (aref vertices-y 1)
                                                             (aref vertices-y 2)
                                                             :bias bias))
                      (intersect-uz-p (triangles-intersect-p (aref vertices-u 0)
                                                             (aref vertices-u 1)
                                                             (aref vertices-u 2)
                                                             (aref vertices-z 0)
                                                             (aref vertices-z 1)
                                                             (aref vertices-z 2)
                                                             :bias bias)))
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
                                                    (make-triangle-mesh (flatten-vertices vertices-w) faces-w
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
                                                    (make-triangle-mesh (flatten-vertices vertices-z) faces-z
                                                                        :diffuse-factor (if intersect-uz-p
                                                                                            #(1 1 0)
                                                                                            #(.7 .7 .6))))
                                              *annotations*)))
                    (org.shirakumo.fraf.wavefront:serialize objects object-file :if-exists :supersede)
                    (render-wavefront object-file output-file :camera-position (dvec 2 2 2))))))))))))

(let ((*material-number* 1)
      (*mesh-number* 1)
      (*annotations* '())
      (*annotation-number* 0)
      (*debug-output* t)

      (bias .125))
  (multiple-value-bind (vertices-u faces-u) (gen-triangle)
    (multiple-value-bind (vertices-v faces-v) (gen-triangle)
      (progn
        (setf (aref vertices-v 0) (aref vertices-u 0)
                                        ; (aref vertices-v 1) (aref vertices-u 1)
              )
        (let* ((u1 (aref vertices-u 0))
               (u2 (aref vertices-u 1))
               (u3 (aref vertices-u 2))
               (v1 (aref vertices-v 0))
               (v2 (aref vertices-v 1))
               (v3 (aref vertices-v 2))
               (intersect-uv-p (triangles-intersect-p u1 u2 u3 v1 v2 v3 :bias bias))
               (v1-on-u-p (point-on-triangle-p u1 u2 u3 v1))
               (v2-on-u-p (point-on-triangle-p u1 u2 u3 v2))
               (v3-on-u-p (point-on-triangle-p u1 u2 u3 v3)))
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
                                            (debug-cube v1 "v1" (make-material (if v1-on-u-p
                                                                                   #(1 0 0)
                                                                                   #(.6 .6 .6)))
                                                        :offset .03)
                                            (debug-cube v2 "v2" (make-material (if v2-on-u-p
                                                                                   #(1 0 0)
                                                                                   #(.6 .6 .6)))
                                                        :offset .03)
                                            (debug-cube v3 "v3" (make-material (if v3-on-u-p
                                                                                   #(1 0 0)
                                                                                   #(.6 .6 .6)))
                                                        :offset .03))
                                      *annotations*)))
            (org.shirakumo.fraf.wavefront:serialize objects object-file :if-exists :supersede)
            (render-wavefront object-file output-file :camera-position (dvec 2 2 2))))))))
