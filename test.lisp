(defpackage #:org.shirakumo.fraf.convex-covering.test
  (:use #:cl #:parachute ; #:org.shirakumo.flare.vector
        )
  (:local-nicknames
   (#:m #:org.shirakumo.fraf.math)
   (#:convex-covering #:org.shirakumo.fraf.convex-covering)
   (#:wavefront #:org.shirakumo.fraf.wavefront))
  (:shadow
   #:test)
  (:export
   #:decompose-file
   #:decompose-test-files))

(in-package #:org.shirakumo.fraf.convex-covering.test)

(defvar *here* #.(make-pathname :name NIL :type NIL :defaults (or *compile-file-pathname* *load-pathname* (error "COMPILE-FILE or LOAD this."))))

(define-test convex-covering)

(defun obj-file (file)
  (merge-pathnames (merge-pathnames file (make-pathname :type "obj" :directory '(:relative "test"))) *here*))

(defun u32 (&rest i)
  (make-array (length i) :element-type '(unsigned-byte 32) :initial-contents i))

(defun f32 (&rest i)
  (make-array (length i) :element-type 'single-float :initial-contents (mapcar #'float i)))

(defun export-hulls (hulls file &key (colors (convex-covering::make-color-generator))
                                     alpha
                                     (highlight nil))
  (wavefront:serialize
   (append (loop for hull across hulls
                 for i :from 0
                 for color = (if colors
                                 (funcall colors i)
                                 (convex-covering::color<-faces
                                  (convex-covering::global-faces hull)))
                 for name = (format NIL "patch~a" i)
                 for mtl = (apply #'make-instance 'wavefront:material
                                  :name name
                                  :diffuse-factor color
                                  (when (and alpha (not (= alpha 1)))
                                    (list :transmission-factor (- 1 alpha))))
                 for mesh = (make-instance 'wavefront:mesh
                                           :name name
                                           :vertex-data (convex-covering:vertices hull)
                                           :index-data (convex-covering:faces hull)
                                           :attributes '(:position)
                                           :material mtl)
                 collect mesh)
           (when highlight
             (let ((bad (make-instance 'wavefront:material
                                       :name "bad"
                                       :diffuse-factor #(0.7 0.2 0.2)))
                   (good (make-instance 'wavefront:material
                                        :name "good"
                                        :diffuse-factor #(0.2 0.7 0.2)))
                   (boundary (make-instance 'wavefront:material
                                            :name "boundary"
                                            :diffuse-factor #(0.7 0.7 0.2)))
                   ;;
                   (patch-edge  (make-instance 'wavefront:material
                                               :name "patch-edge"
                                               :diffuse-factor #(0.7 0.7 0.7)))
                   (facet-edge (make-instance 'wavefront:material
                                              :name "facet-boundary"
                                              :diffuse-factor #(0.7 0.7 0.7)))
                   (facet (make-instance 'wavefront:material
                                         :name "facet"
                                         :diffuse-factor #(0.2 0.2 0.2)))
                   (facet-normal (make-instance 'wavefront:material
                                                :name "facet-normal"
                                                :diffuse-factor #(0.2 0.2 0.7)))
                   ;;
                   (face-normal (make-instance 'wavefront:material
                                               :name "face-normal"
                                               :diffuse-factor #(0.2 0.7 0.7)))
                   (face-edge (make-instance 'wavefront:material
                                             :name "face-edge"
                                             :diffuse-factor #(0.7 0.7 1.0)))
                   (annotation-count 0))
               (loop :with hull = (convex-covering::patch-hull highlight)
                     :for (vertex . kind) :in (remove :face-centroid (convex-covering::annotations hull) :key #'cdr)
                     :for (color offset sample-count)
                        = (multiple-value-list
                           (ecase kind
                             (:bad            (values bad          .1   16))
                             (:good           (values good         .1   16))
                             (:boundary-edge  (values boundary     .01  16))

                             (:patch-edge     (values patch-edge   .01  16))
                             (:facet-edge     (values facet        .01  16))
                             (:facet-normal   (values facet-normal .005 16))
                             (:facet-centroid (values facet-normal .05  16))

                             (:face-normal    (values face-normal  .008 16))
                             (:face-centroid  (values face-normal  .05  16))

                             (:bad-normal     (values bad          .1   16))

                             (:face-edge      (values face-edge    .005 16))))
                     :for name = (format nil "annotation~A" (incf annotation-count))
                     :collect (convex-covering::debug-cube vertex name color :offset offset)))))
   file :if-exists :supersede)
  hulls)

(defun decompose-hulls (file &rest args)
  (let* ((f (wavefront:parse file))
         (meshes (wavefront:extract-meshes f NIL '(:position)))
         (mesh (first meshes))
         (vertices/single-float (wavefront:vertex-data mesh))
         (vertices (map-into (make-array (length vertices/single-float) :element-type 'double-float)
                             (lambda (f) (float f 1.0d0))
                             vertices/single-float))
         (indices (wavefront:index-data mesh)))
    (when (> (length meshes) 1)
      (break "Multiple meshes: ~A" meshes))
    (apply #'convex-covering:decompose vertices indices args)))

(defun decompose-file (in out &rest args)
  (format *trace-output* "; ~A -> ~A~%"
          (enough-namestring in *here*) (enough-namestring out *here*))
  (export-hulls (apply #'decompose-hulls (obj-file in) args)
                (merge-pathnames out (obj-file in))))

(defun decompose-test-files (&rest args)
  (let ((outdir (merge-pathnames (make-pathname :directory '(:relative "test" "output")) *here*)))
    (dolist (in (directory (obj-file (make-pathname :name "gear-1" #+no :wild))))
      (let ((out (make-pathname :name (pathname-name in) :type "obj" :defaults outdir)))
        (time (apply #'decompose-file in out args))))))

(defun test (name &key output (visualization 50))
  (let ((convex-covering::*debug-output* output)
        (convex-covering::*debug-visualizations* visualization))
    (let ((outdir (merge-pathnames (make-pathname :directory '(:relative "test" "output")) *here*)))
      (dolist (in (directory (obj-file (make-pathname :name name #+no :wild))))
        (let ((out (make-pathname :name (pathname-name in) :type "obj" :defaults outdir)))
          (time (decompose-file in out)))))))
