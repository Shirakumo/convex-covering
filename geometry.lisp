(in-package #:org.shirakumo.fraf.convex-covering)

;;; Bounding box

(declaim (ftype (function (dvec3 dvec3) (values vec3 vec3 &optional NIL))
                center-and-size-from-min-and-max)
         ;; (inline center-and-size-from-min-and-max)
         )
(defun center-and-size-from-min-and-max (min max)
  (let* ((size/2 (v* (v- max min) .5))
         (center (v+ min size/2)))
    ;; TODO avoid conversion
    (values (vec (vx center) (vy center) (vz center))
            (vec (vx size/2) (vy size/2) (vz size/2)))))

(declaim (ftype (function (dvec3 dvec3 dvec3) (values vec3 vec3 &optional NIL))
                triangle-bounding-box)
         ;; (inline triangle-bounding-box)
         )
(defun triangle-bounding-box (v1 v2 v3)
  (let ((min (vmin v1 v2 v3))
        (max (vmax v1 v2 v3)))
    (center-and-size-from-min-and-max min max)))

(declaim (ftype (function (manifolds:vertex-array manifolds:face-array manifolds:u32)
                          (values vec3 vec3 &optional NIL))
                face-bounding-box)
         ;; (inline face-bounding-box)
         )
(defun face-bounding-box (vertices faces face-index)
  ;; TODO add manifolds:f
  (let ((v1 (manifolds:v vertices (aref faces (+ (* 3 face-index) 0))))
        (v2 (manifolds:v vertices (aref faces (+ (* 3 face-index) 1))))
        (v3 (manifolds:v vertices (aref faces (+ (* 3 face-index) 2)))))
    (triangle-bounding-box v1 v2 v3)))

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
  (let ((threshold (coerce threshold 'double-float)))
    (cond ((<= (abs (- (min i1-min i2-min) (max i1-max i2-max))) threshold)
           (values T :identical))
          ((<= i1-min i2-min i1-max i2-max) ; case 2
           (d "~4@TOverlap ~A~%"
              (abs (- i2-min i1-max)))
           (values T (if (< (abs (- i2-min i1-max)) threshold) :touching :penetrating)))
          ((<= i1-min i2-min i2-max i1-max) ; case 3
           (d "~4@TI₂ ⊂ I₁ (case 3) Overlaps ~A ~A | threshold ~A~%"
              (abs (- i1-min i2-min))   ; cases 3 4
              (abs (- i1-max i2-max))   ; cases 3 4
              threshold)
           (values T (if (and (< (abs (- i2-min i2-max)) threshold)
                              (or (< (abs (- i1-min i2-min)) threshold)
                                  (< (abs (- i1-max i2-min)) threshold)))
                         :touching
                         :penetrating)))
          ((<= i2-min i1-min i1-max i2-max) ; case 4
           (d "~4@TI₁ ⊂ I₂ (case 4) Overlaps ~A ~A~%"
              (abs (- i1-min i2-min))   ; cases 3 4
              (abs (- i1-max i2-max))   ; cases 3 4
              )
           (values T (if (and (< (abs (- i1-min i1-max)) threshold)
                              (or (< (abs (- i2-min i1-min)) threshold)
                                  (< (abs (- i2-max i1-min)) threshold)))
                         :touching
                         :penetrating)))
          ((<= i2-min i1-min i2-max i1-max) ; case 5
           (d "~4@TOverlap ~A~%"
              (abs (- i1-min i2-max)))
           (values T (if (< (abs (- i1-min i2-max)) threshold) :touching :penetrating)))
          #+no ((not (or (> (- i2-min i1-max) threshold)
                         (> (- i1-min i2-max) threshold)))
                (d "~4@TOverlaps ~A ~A ~A ~A~%"
                   (abs (- i2-min i1-max))  ; case 2
                   (abs (- i1-min i2-min))  ; cases 3 4
                   (abs (- i1-max i2-max))  ; cases 3 4
                   (abs (- i1-min i2-max))) ; case 5
                T)
          (T                            ; no intersection
           NIL))))

(defun triangles-intersect-p (u1 u2 u3 v1 v2 v3 &key (threshold 1d-6))
  (declare (type dvec3 u1 u2 u3 v1 v2 v3))
  (d "U~2@T~A~%~2@T~A~%~2@T~A~%V~2@T~A~%~2@T~A~%~2@T~A~%" u1 u2 u3 v1 v2 v3)
  (flet ((separating-axis-p (axis &optional note)
           (declare (type dvec3 axis))
           (let* ((du1    (v. axis u1))
                  (du2    (v. axis u2))
                  (du3    (v. axis u3))
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
    (multiple-value-bind (result constellation contact)
        (let ((normal-u (vunit (vc (v- u2 u1) (v- u3 u1)))))
          (d "Normal~%")
          (multiple-value-bind (normal-separating-p normal-class)
              (separating-axis-p normal-u 0)
            (d "~2@T=> ~:[not separating~;separating~]~@[ | ~A~]~%"
               normal-separating-p normal-class)
            (cond (normal-separating-p
                   (d "normal separating~%")
                   NIL)
                  ((eq normal-class :identical)
                   (d "Edge normals (coplanar case)~%")
                   (block NIL
                     (let ((contact :penetrating))
                       (flet ((test-axis (e1 e2)
                                (declare (type dvec3 e1 e2))
                                (let ((edge-normal (vunit (vc normal-u (v- e2 e1)))))
                                  (multiple-value-bind (separatingp class)
                                      (separating-axis-p edge-normal 0)
                                    (cond (separatingp
                                           (return (values NIL :coplanar)))
                                          (T
                                           (when (and (eq contact :penetrating)
                                                      (eq class :touching))
                                             (setf contact class))))))))
                         (test-axis u1 u2)
                         (test-axis u2 u3)
                         (test-axis u3 u1)
                         (test-axis v1 v2)
                         (test-axis v2 v3)
                         (test-axis v3 v1)
                         (values T :coplanar contact))))

                   #+old(let ((n12 (vunit (vc normal-u (v- u2 u1))))
                              (n23 (vunit (vc normal-u (v- u3 u2))))
                              (n31 (vunit (vc normal-u (v- u1 u3)))))

                          (cond ((and (not (separating-axis-p n12 0))
                                      (not (separating-axis-p n23 0))
                                      (not (separating-axis-p n31 0)))
                                 (values T :coplanar))
                                (T
                                 (values NIL :coplanar)))))
                  (T
                   (d "General edge pairs~%")
                   (let ((contact normal-class))
                     (if (flet ((test-axis (u1 u2 v1 v2)
                                  (let ((c (vc (v- u2 u1) (v- v2 v1))))
                                    (unless (< (vsqrlength c) threshold)
                                      (let ((c* (vunit c)))
                                        (multiple-value-bind (separatingp class)
                                            (separating-axis-p c* 0)
                                          (when (and (not separatingp)
                                                     (eq contact :penetrating)
                                                     (eq class :touching))
                                            (setf contact class))
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
                                               :do (return-from outer T)))))
                         NIL
                         (values T NIL contact)))))))
      (d "=> ~A~:*~:[~; ~A ~A~]~%"
         result
         (or constellation :general)
         contact)
      (d "--------------------~%")
      (values result constellation contact))))

(defun line-intersects-triangle-p (u1 u2 u3 l1 l2 &key (threshold 1d-6))
  (declare (type dvec3 u1 u2 u3 l1 l2))
  (d "U~2@T~A~%~2@T~A~%~2@T~A~%L~2@T~A~%~2@T~A~%" u1 u2 u3 l1 l2)
  (flet ((separating-axis-p (axis threshold)
           (declare (type dvec3 axis))
           (let* ((du1    (v. axis u1))
                  (du2    (v. axis u2))
                  (du3    (v. axis u3))
                  (dl1    (v. axis l1))
                  (dl2    (v. axis l2))
                  (iu-min (min du1 du2 du3)) ; interval u min
                  (iu-max (max du1 du2 du3)) ; interval u max
                  (iv-min (min dl1 dl2))     ; interval v min
                  (iv-max (max dl1 dl2)))
             (d "~2@TAxis~@[ ~A~] ~A~&~
                 ~4@T=> {~5,2F ~5,2F ~5,2F} {~5,2F ~5,2F}~&~
                 ~4@T=> [~5,2F, ~5,2F] [~5,2F, ~5,2F]~%"
                NIL axis
                du1 du2 du3 dl1 dl2
                iu-min iu-max iv-min iv-max)
             (multiple-value-bind (intersectp class)
                 (intervals-intersect-p iu-min iu-max iv-min iv-max
                                        :threshold threshold)
               (d "~4@T=> ~A~@[ | ~A~]~%" (if intersectp "not separating" "separating") class)
               (values (not intersectp) class)))))
    (d "--------------------~%")
    (multiple-value-bind (result constellation contact)
        (let ((normal-u (vunit (vc (v- u2 u1) (v- u3 u1)))))
          (d "Normal~%")
          (multiple-value-bind (normal-separating-p normal-class)
              (separating-axis-p normal-u 1d-3) ; TODO
            (d "~2@T=> ~:[not separating~;separating~]~@[ | ~A~]~%"
               normal-separating-p normal-class)
            (cond (normal-separating-p
                   (d "normal separating~%")
                   NIL)
                  ((eq normal-class :identical)
                   (d "Edge normals (coplanar case)~%")
                   (block NIL
                     (let ((contact :penetrating))
                       (flet ((test-axis (e1 e2)
                                (declare (type dvec3 e1 e2))
                                (let ((edge-normal (vunit (vc normal-u (v- e2 e1)))))
                                  (multiple-value-bind (separatingp class)
                                      (separating-axis-p edge-normal threshold)
                                    (cond (separatingp
                                           (return (values NIL :coplanar)))
                                          ((and (eq class :touching)
                                                (eq contact :penetrating))
                                           (setf contact class)))))))
                         (test-axis u1 u2)
                         (test-axis u2 u3)
                         (test-axis u3 u1)
                         (test-axis l1 l2)
                         (if (and (eq contact :penetrating)
                                  (member l1 (list u1 u2 u3) :test (lambda (v1 v2)
                                                                     (< (vsqrdistance v1 v2) (* threshold threshold))))
                                  (member l2 (list u1 u2 u3) :test (lambda (v1 v2)
                                                                     (< (vsqrdistance v1 v2) (* threshold threshold)))))
                             (values T :coplanar :edge)
                             (values T :coplanar contact))))))
                  (T
                   (d "General edge pairs~%")
                   (let ((contact normal-class))
                     (if (flet ((test-axis (u1 u2 v1 v2)
                                  (let ((c (vc (v- u2 u1) (v- v2 v1))))
                                    (unless (< (vsqrlength c) (* threshold threshold))
                                      (let ((c* (vunit c)))
                                        (multiple-value-bind (separatingp class)
                                            (separating-axis-p c* threshold)
                                          (when (and (not separatingp)
                                                     (eq contact :penetrating)
                                                     (eq class :touching))
                                            (setf contact class))
                                          (when (boundp '*annotation-number*)
                                            (apply #'debug-line** u1 (v* c* .1)
                                                   (when separatingp
                                                     (list :diffuse-factor #(1 0 0)))))
                                          separatingp))))))
                           ;; TODO(jmoringe): The last call is a hack
                           ;; in case the triangle and the line
                           ;; segment are actually coplanar but we
                           ;; didn't detect it due to insufficient
                           ;; precision.
                           (or (loop for (a1 a2) on (list u1 u2 u3 u1)
                                     while a2
                                     do (when (test-axis a1 a2 l1 l2)
                                          (return T)))
                               (and (> (vsqrlength (vc normal-u (v- l1 l2))) (* threshold threshold))
                                    (separating-axis-p (vunit (vc normal-u (v- l1 l2))) threshold))))
                         NIL
                         (values T NIL contact)))))))
      (d "=> ~A~:*~:[~; ~A ~A~]~%"
         result
         (or constellation :general)
         contact)
      (d "--------------------~%")
      (values result constellation contact))))

(defun point-on-triangle-p (u1 u2 u3 v &key (threshold .00001))
  (let ((contact :separating))
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
               (multiple-value-bind (result contact*)
                   (intervals-intersect-p iu-min iu-max dv dv
                                          :threshold threshold)
                 (when (eq contact* :touching)
                   (setf contact contact*))
                 (d "~2@T~A~%" (if result "not separating" "separating"))
                 (not result)))))
      (let* ((normal-u (vunit (vc (v- u2 u1) (v- u3 u1))))
             (n12      (vunit (vc normal-u (v- u2 u1)))) ; TODO compute when needed
             (n23      (vunit (vc normal-u (v- u3 u2))))
             (n31      (vunit (vc normal-u (v- u1 u3)))))
        (d "--------------------~%")
        (multiple-value-prog1
            (cond ((not (< (abs (v. normal-u (v- v u1))) threshold))
                   nil)
                  ((and (not (separating-axis-p n12 0))
                        (not (separating-axis-p n23 0))
                        (not (separating-axis-p n31 0)))
                   (values t contact)))
          (d "--------------------~%"))))))
