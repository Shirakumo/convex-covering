(defpackage #:org.shirakumo.fraf.convex-covering.test
  (:use #:cl #:parachute #:org.shirakumo.flare.vector)
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

(defun export-hulls (hulls file)
  (wavefront:serialize
   (loop with colors = #(#(1.0 0.0 0.0)
                         #(0.0 1.0 0.0)
                         #(0.0 0.0 1.0)
                         #(0.5 0.5 0.0)
                         #(0.0 0.5 0.5)
                         #(0.5 0.0 0.5)
                         #(1.0 1.0 1.0)
                         #(0.1 0.1 0.1))
         for hull across hulls
         for i from 0
         for mtl = (make-instance 'wavefront:material
                                  :name (format NIL "c~a" i)
                                  :diffuse-factor (aref colors (mod i (length colors))))
         collect (make-instance 'wavefront:mesh
                                :vertex-data (convex-covering:vertices hull)
                                :index-data (convex-covering:faces hull)
                                :attributes '(:position)
                                :material mtl))
   file :if-exists :supersede)
  hulls)

(defun decompose-hulls (file &rest args)
  (let* ((f (wavefront:parse file))
         (m (first (wavefront:extract-meshes f NIL '(:position)))))
    (apply #'convex-covering::decompose
           (wavefront:vertex-data m)
           (wavefront:index-data m)
           args)))

(defun decompose-file (in out &rest args)
  (export-hulls (apply #'decompose-hulls (obj-file in) args)
                (merge-pathnames out (obj-file in))))

(defun decompose-test-files (&rest args)
  (let ((outdir (merge-pathnames (make-pathname :directory '(:relative "test" "output")) *here*)))
    (dolist (in (directory (test-file (make-pathname :name :wild))))
      (let ((out (make-pathname :name (pathname-name in) :type "obj" :defaults outdir)))
        (apply #'decompose-file in out args)))))
