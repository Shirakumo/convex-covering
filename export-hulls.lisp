(in-package #:org.shirakumo.fraf.convex-covering)

(defun export-hulls (hulls file &key (colors (make-color-generator))
                                     alpha
                                     (highlight nil))
  (org.shirakumo.fraf.wavefront:serialize
   (append (loop for hull across hulls
                 for i :from 0
                 for color = (if colors
                                 (funcall colors i)
                                 (color<-faces (etypecase hull
                                                 (hull (hull-global-faces hull))
                                                 (convex-hull (faces hull)))))
                 for name = (format NIL "patch~a" i)
                 for mtl = (apply #'make-instance 'org.shirakumo.fraf.wavefront:material
                                  :name name
                                  :diffuse-factor color
                                  (when (and alpha (not (= alpha 1)))
                                    (list :transmission-factor (- 1 alpha))))
                 for mesh = (make-instance 'org.shirakumo.fraf.wavefront:mesh
                                           :name name
                                           :vertex-data (etypecase hull
                                                          (hull (hull-vertices hull))
                                                          (convex-hull (vertices hull)))
                                           :index-data (etypecase hull
                                                         (hull (hull-facets hull))
                                                         (convex-hull (faces hull)))
                                           :attributes '(:position)
                                           :material mtl)
                 collect mesh)
           (when highlight
             (let ((bad (make-instance 'org.shirakumo.fraf.wavefront:material
                                       :name "bad"
                                       :diffuse-factor #(0.7 0.2 0.2)))
                   (good (make-instance 'org.shirakumo.fraf.wavefront:material
                                        :name "good"
                                        :diffuse-factor #(0.2 0.7 0.2)))
                   (boundary (make-instance 'org.shirakumo.fraf.wavefront:material
                                            :name "boundary"
                                            :diffuse-factor #(0.7 0.7 0.2)))
                   ;;
                   (patch-edge  (make-instance 'org.shirakumo.fraf.wavefront:material
                                               :name "patch-edge"
                                               :diffuse-factor #(0.7 0.7 0.7)))
                   (facet-edge (make-instance 'org.shirakumo.fraf.wavefront:material
                                              :name "facet-boundary"
                                              :diffuse-factor #(0.7 0.7 0.7)))
                   (facet (make-instance 'org.shirakumo.fraf.wavefront:material
                                         :name "facet"
                                         :diffuse-factor #(0.2 0.2 0.2)))
                   (facet-normal (make-instance 'org.shirakumo.fraf.wavefront:material
                                                :name "facet-normal"
                                                :diffuse-factor #(0.2 0.2 0.7)))
                   ;;
                   (face-normal (make-instance 'org.shirakumo.fraf.wavefront:material
                                               :name "face-normal"
                                               :diffuse-factor #(0.2 0.7 0.7)))
                   (face-edge (make-instance 'org.shirakumo.fraf.wavefront:material
                                             :name "face-edge"
                                             :diffuse-factor #(0.7 0.7 1.0)))
                   (annotation-count 0))
               (loop :with hull = (patch-hull highlight)
                     :for (vertex . kind) :in (remove :face-centroid (hull-annotations hull) :key #'cdr)
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
                     :collect (debug-cube vertex name color :offset offset)))))
   file :if-exists :supersede)
  hulls)
