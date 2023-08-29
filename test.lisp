(defpackage #:org.shirakumo.fraf.convex-covering.test
  (:use #:cl #:parachute ; #:org.shirakumo.flare.vector
        )
  (:local-nicknames
   (#:convex-covering #:org.shirakumo.fraf.convex-covering)
   (#:wavefront #:org.shirakumo.fraf.wavefront))
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

(defun export-hulls (hulls file &key (colors (make-color-generator)))
  (wavefront:serialize
   (loop for hull across hulls
         for i :from 0
         for color = (funcall colors i)
         for name = (format NIL "patch~a" i)
         for mtl = (make-instance 'wavefront:material
                                  :name name
                                  :diffuse-factor color)
         for mesh = (make-instance 'wavefront:mesh
                                   :name name
                                   :vertex-data (convex-covering:vertices hull)
                                   :index-data (convex-covering:faces hull)
                                   :attributes '(:position)
                                   :material mtl)
         collect mesh)
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
  (export-hulls (apply #'decompose-hulls (obj-file in) args)
                (merge-pathnames out (obj-file in))))

(defun decompose-test-files (&rest args)
  (let ((outdir (merge-pathnames (make-pathname :directory '(:relative "test" "output")) *here*)))
    (dolist (in (directory (obj-file (make-pathname :name "gear-2" #+no :wild))))
      (let ((out (make-pathname :name (pathname-name in) :type "obj" :defaults outdir)))
        (format *trace-output* "; ~A -> ~A~%"
                (enough-namestring in *here*) (enough-namestring out *here*))
        (time (apply #'decompose-file in out args))))))
