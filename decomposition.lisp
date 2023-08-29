;;;; This is an implementation of the algorithm from the paper
;;;;   Convex Hull Covering of Polygonal Scenes for Accurate Collision Detection in Games
;;;; by Rong Liu et al. accessible at https://www.cs.sfu.ca/~haoz/pubs/liu_zhang_gi08.pdf
;;;;

(in-package #:org.shirakumo.fraf.convex-covering)

;; TODO cache face normals in this structure?
(defstruct (convex-hull
            (:constructor make-convex-hull (vertices faces))
            (:conc-name NIL)
            (:copier NIL)
            (:predicate NIL))
  (vertices  #() :type (manifolds:vertex-array double-float))
  (faces #() :type manifolds:face-array))

(defstruct (patch
            (:constructor %make-patch (faces &optional hull compactness links))
            (:copier NIL)
            (:predicate NIL))
  (faces #() :type manifolds:face-array)
  (hull NIL :type (or null convex-hull))
  (links (make-array 0 :adjustable T :fill-pointer T) :type vector)
  (compactness 0.0d0 :type double-float))

(defun make-patch (a b c)
  (%make-patch (make-array 3 :element-type '(unsigned-byte 32) :initial-contents (list a b c))))

(defun compute-patch-convex-hull (all-vertices faces)
  (let* ((vertex-count (length faces))
         (vertices ; (make-array (* 3 vertex-count) :element-type 'double-float)
           (make-array 9 :element-type 'double-float :adjustable t :fill-pointer 0)))
    (let ((seen (make-hash-table :test #'equal)))
      (loop :for i :below vertex-count
            :for j = (aref faces i)
            :for key = (list (aref all-vertices (+ (* 3 j) 0))
                             (aref all-vertices (+ (* 3 j) 1))
                             (aref all-vertices (+ (* 3 j) 2)))
            :do (unless (gethash key seen) ; in case quickhull does not permit duplicate vertices
                  (setf (gethash key seen) t)
                  #+no (setf (aref vertices (+ (* 3 i) 0)) (aref all-vertices (+ (* 3 j) 0))
                             (aref vertices (+ (* 3 i) 1)) (aref all-vertices (+ (* 3 j) 1))
                             (aref vertices (+ (* 3 i) 2)) (aref all-vertices (+ (* 3 j) 2)))
                  (vector-push-extend (aref all-vertices (+ (* 3 j) 0)) vertices)
                  (vector-push-extend (aref all-vertices (+ (* 3 j) 1)) vertices)
                  (vector-push-extend (aref all-vertices (+ (* 3 j) 2)) vertices))))
    #+no (map-into vertices (lambda (index)
                              (aref all-vertices index)))
    (multiple-value-call #'make-convex-hull
      (org.shirakumo.fraf.quickhull:convex-hull vertices))))

(defun push-link (new-link patch)
  (assert (or (eq patch (patch-link-a new-link))
              (eq patch (patch-link-b new-link))))
  (vector-push-extend new-link (patch-links patch)))

(defun find-link-involving (other-patch patch)
  (find-if (lambda (link)
             (link-involves-p other-patch link))
           (patch-links patch)))

;;;

(defun make-merged-patch (all-vertices all-faces a b)
  (let* ((faces1 (patch-faces a))
         (faces2 (patch-faces b))
         (faces (make-array (+ (length faces1) (length faces2)) :element-type '(unsigned-byte 32)))) ; TODO just `concatenate'?
    (replace faces faces1)
    (replace faces faces2 :start1 (length faces1))
    (let* ((hull (compute-patch-convex-hull all-vertices faces))
           (patch (%make-patch faces hull)))
      (when (valid-patch-p patch all-vertices all-faces)
        (setf (patch-compactness patch) (compute-compactness all-vertices patch))
        patch))))

(defstruct (patch-link
            (:constructor %make-patch-link (a b &optional merge-cost merge-result))
            (:copier NIL)
            (:predicate NIL))
  (a NIL :type patch)
  (b NIL :type patch)
  (merge-cost most-positive-double-float :type double-float)
  (merge-result NIL :type (or null patch)))

(defun make-patch-link (all-vertices all-faces patch1 patch2)
  (assert (not (eq patch1 patch2)))
  (let ((result (make-merged-patch all-vertices all-faces patch1 patch2)))
    (if result
        (%make-patch-link patch1 patch2 (/ (patch-compactness result)) result)
        (%make-patch-link patch1 patch2))))

(defun link-patches (all-vertices all-faces a b) ; called in preparation phase for adjacent faces
  (unless (loop for link across (patch-links a) ; TODO make a predicate
                thereis (or (eql b (patch-link-a link))
                            (eql b (patch-link-b link))))
    (let ((link (make-patch-link all-vertices all-faces a b)))
      (push-link link a)
      (push-link link b)
      link)))

(defun link-involves-p (patch link)
  (or (eq (patch-link-a link) patch)
      (eq (patch-link-b link) patch)))

(defun merge-patches (all-vertices all-faces link)
  (let ((patch (patch-link-merge-result link)))
    ;; Now that we actually merge this in, compute new links and update the existing neighbour's links.
    (flet ((link (cur other)
             (loop with links = (alexandria:copy-array (patch-links cur))
                   for i from 0 below (length links)
                   for link = (aref links i)
                   do (unless (link-involves-p other link) #+was (or (eql other (patch-link-a link))
                                                                     (eql other (patch-link-b link)))
                        (let* ((patch2 (link-other-patch link cur) #+no (if (eql cur (patch-link-a link))
                                                                       (patch-link-b link)
                                                                       (patch-link-a link)))
                               (new-link (make-patch-link all-vertices all-faces patch patch2)))
                          (unless (find-link-involving patch2 patch)
                            (push-link new-link patch))
                          (let ((temp (remove-if (lambda (old-link)
                                                   (or (eq old-link link)
                                                       (link-involves-p cur old-link)
                                                       (link-involves-p patch old-link)))
                                                 (patch-links patch2))))
                           (setf (patch-links patch2)
                                 (make-array (length temp) :initial-contents temp :adjustable t :fill-pointer t)))
                          (push-link new-link patch2)
                          #+no (let ((j (position link (patch-links patch2))))
                                 (setf (aref (patch-links patch2) j) new-link)))))))
      (link (patch-link-a link) (patch-link-b link))
      (link (patch-link-b link) (patch-link-a link))
      patch)))

(defvar *winner* nil)

(defun patch-debug-name (patch)
  (let* ((faces (patch-faces patch))
         (count (length faces)))
    (format nil "~{[~D ~D ~D]~^ ~}~@[ …~]"
            (coerce (subseq faces 0 (min 9 count)) 'list)
            (> count 9))))

(defun next-link (links)
  (loop with min = most-positive-double-float
        with min-link = NIL
        for link being the hash-keys of links
        for cost = (patch-link-merge-cost link)
        #+no :do #+no (format *trace-output* "?? ~5,2F ~{[~D ~D ~D]~^ ~}/~A~%  -- ~{[~D ~D ~D]~^ ~}/~A~%"
                    cost
                    (coerce (patch-faces (patch-link-a link)) 'list)
                    (manifolds:boundary-list (patch-faces (patch-link-a link)))
                    (coerce (patch-faces (patch-link-b link)) 'list)
                    (manifolds:boundary-list (patch-faces (patch-link-b link))))
        do (when (<= cost min)
             (setf min (patch-link-merge-cost link))
             (setf min-link link))
        finally (when min-link
                  (format *trace-output* "=> ~5,2F ~A -- ~A~%"
                          (patch-link-merge-cost min-link)
                          (patch-debug-name (patch-link-a min-link))
                          (patch-debug-name (patch-link-b min-link))))
                (setf *winner* min-link)
                (return min-link)))

(defun link-other-patch (link this-patch)
  (cond ((eq this-patch (patch-link-a link))
         (patch-link-b link))
        ((eq this-patch (patch-link-b link))
         (patch-link-a link))
        (t
         (error "corrupt link"))))

(defun decompose (vertices indices &key) ; TODO indices -> faces
  ;; FIXME: This is all really dumb and uses really bad data structures
  ;;        Could definitely be optimised a lot by someone smarter
  (let ((patchlist (make-array (truncate (length indices) 3)))
        (patches (make-hash-table :test 'eq))
        (links (make-hash-table :test 'eq))
        (i 0))
    ;; 1. Destructure the mesh into one patch per face
    (manifolds:do-faces (a b c indices)
      (let ((patch (make-patch a b c)))
                                        ; (setf (patch-hull patch) (compute-patch-convex-hull vertices (patch-faces patch)))
        (setf (gethash patch patches) T)
        (setf (aref patchlist i) patch)
        (incf i)))
    ;; 2. Find neighbouring patches and create the links
    (let ((adjacents (manifolds:face-adjacency-list indices)))
      (dotimes (face (length adjacents))
        (loop for other in (aref adjacents face)
              for link = (link-patches vertices indices (aref patchlist face) (aref patchlist other))
              when link
              do (setf (gethash link links) T))))
    ;; 3. Greedily merge patches according to merge cost
    (loop                               ; :repeat 7
          :for i :from 1
          for link = (next-link links)
          while (and link (patch-link-merge-result link)) ; TODO for after while is not conforming
          for patch1 = (patch-link-a link)
          for patch2 = (patch-link-b link)

          :do (format *trace-output* "------------ Step ~D ------------~%" i)
              (let ((patches (alexandria:hash-table-keys patches)))
                (visualize-step patches i))

          do ;; 1. Remove old patches and links
             (format *trace-output* "~:D patch~:P ~D link~:P~%"
                     (hash-table-count patches) (hash-table-count links))
             (remhash patch1 patches)
             (remhash patch2 patches)
             (remhash link links)
             (flet ((remove-other-links (patch)
                      (loop for other-link across (patch-links patch)
                                        ; for other-patch = (link-other-patch other-link patch)
                            do (remhash other-link links)
                                        ; (setf (patch-links other-patch) (remove other-link (patch-links other-patch)))
                            )))
               (remove-other-links patch1)
               (remove-other-links patch2))
             (format *trace-output* "  after removing ~:D patch~:P ~D link~:P~%"
                     (hash-table-count patches) (hash-table-count links))
             ;; 2. Merge the patches
             (let ((patch (merge-patches vertices indices link)))
               ;; 3. Insert the new links
               (assert (patch-hull patch))
               (setf (gethash patch patches) T)
               (loop for link across (patch-links patch)
                     do (setf (gethash link links) T))))
    ;; 4. Return the patches' convex hulls
    (let ((hulls (make-array (hash-table-count patches))))
      (loop for patch being the hash-keys of patches
            for i from 0
            do (setf (aref hulls i) (patch-hull patch)))
      (remove nil hulls)                ; TODO temp hack
      )))

;;; Utilities

(defun debug-filename (prefix i type)
  (format nil "/tmp/~A-~3,'0D.~A" prefix i type))

;;; Graph

(defclass decomposition ()
  ((%colors :initarg :colors
            :reader  colors)
   (%order  :initarg :order
            :reader  order)))

(defmethod cl-dot:graph-object-node ((graph decomposition) (object patch))
  (let* ((index (position object (order graph)))
         (color (when index
                  (funcall (colors graph) index)))
         (label (patch-debug-name object)))
    (make-instance 'cl-dot:node :attributes `(:label ,label
                                              ,@(when color
                                                  `(:color ,(format nil "#~2,'0X~2,'0X~2,'0X"
                                                                    (floor (aref color 0) 1/255)
                                                                    (floor (aref color 1) 1/255)
                                                                    (floor (aref color 2) 1/255))))))))

(defmethod cl-dot:graph-object-points-to ((graph decomposition) (object patch))
  (flet ((make-edge (link)
           (let ((color (if (patch-link-merge-result link)
                            "green"
                            "red"))
                 (width (if (eq link *winner*) 5 nil)))
             (make-instance 'cl-dot:attributed :object     (link-other-patch link object)
                                               :attributes `(:color ,color
                                                             ,@(when width `(:penwidth ,width)))))))
    (map 'list #'make-edge (patch-links object))))

(defun graph-step (patches i &key (output-file (debug-filename "graph" i "png"))
                                  (colors      (org.shirakumo.fraf.convex-covering.test::make-color-generator)))
  (let* ((client (make-instance 'decomposition :colors colors :order (remove nil patches :key #'patch-hull)))
         (graph  (cl-dot:generate-graph-from-roots client patches)))
    (cl-dot:dot-graph graph output-file :format :png)))

;;; Render

(defun render-step (patches i &key (output-file (debug-filename "hulls" i "png"))
                                   (colors      (org.shirakumo.fraf.convex-covering.test::make-color-generator)))
  (let ((hulls           (remove nil (map 'vector #'patch-hull patches)))
        (object-filename (debug-filename "hulls" i "obj")))
    (org.shirakumo.fraf.convex-covering.test::export-hulls
     hulls object-filename :colors colors)
    (inferior-shell:run `("f3d" "--output" ,output-file
                                "--camera-position" "5,4,5"
                                ,object-filename))))

;;; Visualization

(defun visualize-step (patches i)
  (let ((graph-filename (debug-filename "graph" i "png"))
        (image-filename (debug-filename "hulls" i "png"))
        (step-filename  (debug-filename "step" i "png")))
    ; (graph-step patches i :output-file graph-filename)
    (render-step patches i :output-file image-filename)
    #+no (inferior-shell:run `("montage" ,graph-filename ,image-filename
                                    "-tile" "2x1" "-geometry" "+0+0"
                                    ,step-filename))))
