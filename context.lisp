(in-package #:org.shirakumo.fraf.convex-covering)

;;; Vertex index
;;;
;;; Lookup index in global vertex array given the vertex coordinates.

(defun make-vertex-index ()
  (make-hash-table :test #'equalp))

(defun vertex-position (vertex index)
  (check-type vertex dvec3)
  (the (or null manifolds:u32) (gethash vertex index)))

(defun (setf vertex-position) (new-value vertex index)
  (check-type vertex dvec3)
  (setf (gethash vertex index) new-value))

(defun index-vertices (vertices)
  (check-type vertices manifolds:vertex-array)
  (loop with index = (make-vertex-index)
        for i below (/ (length vertices) 3)
        for vertex = (manifolds:v vertices i)
        do (setf (vertex-position vertex index) i)
        finally (return index)))

;;; Face index

(defstruct (face-info
            (:constructor make-face-info (index center size/2 normal)))
  (index  (error "required") :type manifolds:u32 :read-only T)
  (center (error "required") :type vec3          :read-only T)
  (size/2 (error "required") :type vec3          :read-only T)
  (normal (error "required") :type dvec3         :read-only T))

(defmethod space:location ((object face-info))
  (face-info-center object))

(defmethod space:bsize ((object face-info))
  (face-info-size/2 object))

(defun index-faces (vertices faces)
  (let ((index (ecase :kd-tree
                 (:grid    (org.shirakumo.fraf.trial.space.grid3:make-grid .3 :bsize (vec3 2 2 2)))
                 (:kd-tree (org.shirakumo.fraf.trial.space.kd-tree:make-kd-tree :dimensions 3)))))
    (loop for i below (/ (length faces) 3)
          for (center size/2) = (multiple-value-list ; SBCL optimizes this away
                                 (face-bounding-box vertices faces i))
          for normal = (vunit (manifolds:face-normal vertices faces i))
          for info = (make-face-info i center size/2 normal)
          do (space:enter info index))
    (space:reoptimize index) ; TODO use enter-all
    index))

;;; Boundary edge index

(defstruct (boundary-bar-info
            (:constructor make-boundary-bar-info (bar1 bar2)))
  (bar1 (error "required") :type dvec3 :read-only t)
  (bar2 (error "required") :type dvec3 :read-only t))

(defmethod space:location ((object boundary-bar-info))
  (let ((a (boundary-bar-info-bar1 object))
        (b (boundary-bar-info-bar2 object)))
    ;; TODO avoid conversion
    (vec (v/ (v+ a b) 2))))

(defmethod space:bsize ((object boundary-bar-info))
  (let ((a (boundary-bar-info-bar1 object))
        (b (boundary-bar-info-bar2 object)))
    ;; TODO avoid conversion
    (vec (v/ (vabs (v- a b)) 2))))

(defun compute-boundary-bar (all-vertices boundary-edge)
  (let* ((a (manifolds:start boundary-edge))
         (b (manifolds:end boundary-edge))
         (c (manifolds:opposite boundary-edge))
         (v1 (manifolds:v all-vertices a))
         (v2 (manifolds:v all-vertices b))
         (v3 (manifolds:v all-vertices c))
         (edge-center (v/ (v+ v1 v2) 2))
         (face-normal (vunit (vc (v- v2 v1) (v- v3 v1))))
         ;; TODO(jmoringe): the paper doesn't say how to scale this
         (offset (v* face-normal .1))
         (bar1 (v- edge-center offset))
         (bar2 (v+ edge-center offset)))
    (values bar1 bar2)))

(defun index-boundary-edges (all-vertices boundary-edges)
  (let ((index (org.shirakumo.fraf.trial.space.kd-tree:make-kd-tree :dimensions 3)))
    (loop for boundary-edge :across boundary-edges
          for (bar1 bar2) = (multiple-value-list (compute-boundary-bar
                                                  all-vertices boundary-edge))
          for info = (make-boundary-bar-info bar1 bar2)
          do (space:enter info index))
    (space:reoptimize index) ; TODO enter-all
    index))

;;; Context

(defstruct (context
            (:constructor make-context (vertices faces
                                        &aux (boundary-edges      (manifolds::boundary-list faces))
                                             (vertex-index        (index-vertices vertices))
                                             (face-index          (index-faces vertices faces))
                                             (boundary-edge-index (index-boundary-edges
                                                                   vertices boundary-edges))))
            (:copier nil)
            (:predicate nil))
  ;; Mesh
  (vertices            (error "required") :type (manifolds:vertex-array double-float) :read-only T)
  (faces               (error "required") :type manifolds:face-array :read-only T)
  (boundary-edges      (error "required") :read-only T) ; TODO types
  ;; Index structures
  (vertex-index        (error "required") :read-only T)
  (face-index          (error "required") :read-only T)
  (boundary-edge-index (error "required") :read-only T))